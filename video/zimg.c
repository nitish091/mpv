/*
 * This file is part of mpv.
 *
 * mpv is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * mpv is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with mpv.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <math.h>

#include <libavutil/bswap.h>
#include <libavutil/pixfmt.h>

#include "common/common.h"
#include "common/msg.h"
#include "csputils.h"
#include "options/m_config.h"
#include "options/m_option.h"
#include "repack.h"
#include "video/fmt-conversion.h"
#include "video/img_format.h"
#include "zimg.h"

static_assert(MP_IMAGE_BYTE_ALIGN >= ZIMG_ALIGN, "");

#define HAVE_ZIMG_ALPHA (ZIMG_API_VERSION >= ZIMG_MAKE_API_VERSION(2, 4))

static const struct m_opt_choice_alternatives mp_zimg_scalers[] = {
    {"point",           ZIMG_RESIZE_POINT},
    {"bilinear",        ZIMG_RESIZE_BILINEAR},
    {"bicubic",         ZIMG_RESIZE_BICUBIC},
    {"spline16",        ZIMG_RESIZE_SPLINE16},
    {"spline36",        ZIMG_RESIZE_SPLINE36},
    {"lanczos",         ZIMG_RESIZE_LANCZOS},
    {0}
};

#define OPT_PARAM(var) OPT_DOUBLE(var), .flags = M_OPT_DEFAULT_NAN

#define OPT_BASE_STRUCT struct zimg_opts
const struct m_sub_options zimg_conf = {
    .opts = (struct m_option[]) {
        {"scaler", OPT_CHOICE_C(scaler, mp_zimg_scalers)},
        {"scaler-param-a", OPT_PARAM(scaler_params[0])},
        {"scaler-param-b", OPT_PARAM(scaler_params[1])},
        {"scaler-chroma", OPT_CHOICE_C(scaler_chroma, mp_zimg_scalers)},
        {"scaler-chroma-param-a", OPT_PARAM(scaler_chroma_params[0])},
        {"scaler-chroma-param-b", OPT_PARAM(scaler_chroma_params[1])},
        {"dither", OPT_CHOICE(dither,
            {"no",              ZIMG_DITHER_NONE},
            {"ordered",         ZIMG_DITHER_ORDERED},
            {"random",          ZIMG_DITHER_RANDOM},
            {"error-diffusion", ZIMG_DITHER_ERROR_DIFFUSION})},
        {"fast", OPT_FLAG(fast)},
        {0}
    },
    .size = sizeof(struct zimg_opts),
    .defaults = &(const struct zimg_opts){
        .scaler = ZIMG_RESIZE_LANCZOS,
        .scaler_params = {NAN, NAN},
        .scaler_chroma_params = {NAN, NAN},
        .scaler_chroma = ZIMG_RESIZE_BILINEAR,
        .dither = ZIMG_DITHER_RANDOM,
        .fast = 1,
    },
};

struct mp_zimg_repack {
    bool pack;                  // if false, this is for unpacking
    struct mp_image_params fmt; // original mp format (possibly packed format,
                                // swapped endian)
    int zimgfmt;                // zimg equivalent unpacked format
    int num_planes;             // number of planes involved
    unsigned zmask[4];          // zmask[mp_index] = zimg mask (using mp index!)
    int z_planes[4];            // z_planes[zimg_index] = mp_index (or -1)

    struct mp_repack *repack;   // converting to/from planar

    // Temporary memory for slice-wise repacking. This may be set even if repack
    // is not set (then it may be used to avoid alignment issues). This has
    // about one slice worth of data.
    struct mp_image *tmp;

    int real_w, real_h;         // aligned size
};

static void mp_zimg_update_from_cmdline(struct mp_zimg_context *ctx)
{
    m_config_cache_update(ctx->opts_cache);

    struct zimg_opts *opts = ctx->opts_cache->opts;
    ctx->opts = *opts;
}

static zimg_chroma_location_e mp_to_z_chroma(enum mp_chroma_location cl)
{
    switch (cl) {
    case MP_CHROMA_LEFT:        return ZIMG_CHROMA_LEFT;
    case MP_CHROMA_CENTER:      return ZIMG_CHROMA_CENTER;
    default:                    return ZIMG_CHROMA_LEFT;
    }
}

static zimg_matrix_coefficients_e mp_to_z_matrix(enum mp_csp csp)
{
    switch (csp) {
    case MP_CSP_BT_601:         return ZIMG_MATRIX_BT470_BG;
    case MP_CSP_BT_709:         return ZIMG_MATRIX_BT709;
    case MP_CSP_SMPTE_240M:     return ZIMG_MATRIX_ST240_M;
    case MP_CSP_BT_2020_NC:     return ZIMG_MATRIX_BT2020_NCL;
    case MP_CSP_BT_2020_C:      return ZIMG_MATRIX_BT2020_CL;
    case MP_CSP_RGB:            return ZIMG_MATRIX_RGB;
    case MP_CSP_XYZ:            return ZIMG_MATRIX_RGB;
    case MP_CSP_YCGCO:          return ZIMG_MATRIX_YCGCO;
    default:                    return ZIMG_MATRIX_BT709;
    }
}

static zimg_transfer_characteristics_e mp_to_z_trc(enum mp_csp_trc trc)
{
    switch (trc) {
    case MP_CSP_TRC_BT_1886:    return ZIMG_TRANSFER_BT709;
    case MP_CSP_TRC_SRGB:       return ZIMG_TRANSFER_IEC_61966_2_1;
    case MP_CSP_TRC_LINEAR:     return ZIMG_TRANSFER_LINEAR;
    case MP_CSP_TRC_GAMMA22:    return ZIMG_TRANSFER_BT470_M;
    case MP_CSP_TRC_GAMMA28:    return ZIMG_TRANSFER_BT470_BG;
    case MP_CSP_TRC_PQ:         return ZIMG_TRANSFER_ST2084;
    case MP_CSP_TRC_HLG:        return ZIMG_TRANSFER_ARIB_B67;
    case MP_CSP_TRC_GAMMA18:    // ?
    case MP_CSP_TRC_GAMMA20:
    case MP_CSP_TRC_GAMMA24:
    case MP_CSP_TRC_GAMMA26:
    case MP_CSP_TRC_PRO_PHOTO:
    case MP_CSP_TRC_V_LOG:
    case MP_CSP_TRC_S_LOG1:
    case MP_CSP_TRC_S_LOG2:     // ?
    default:                    return ZIMG_TRANSFER_BT709;
    }
}

static zimg_color_primaries_e mp_to_z_prim(enum mp_csp_prim prim)
{
    switch (prim) {
    case MP_CSP_PRIM_BT_601_525:return ZIMG_PRIMARIES_ST170_M;
    case MP_CSP_PRIM_BT_601_625:return ZIMG_PRIMARIES_BT470_BG;
    case MP_CSP_PRIM_BT_709:    return ZIMG_PRIMARIES_BT709;
    case MP_CSP_PRIM_BT_2020:   return ZIMG_PRIMARIES_BT2020;
    case MP_CSP_PRIM_BT_470M:   return ZIMG_PRIMARIES_BT470_M;
    case MP_CSP_PRIM_CIE_1931:  return ZIMG_PRIMARIES_ST428;
    case MP_CSP_PRIM_DCI_P3:    return ZIMG_PRIMARIES_ST431_2;
    case MP_CSP_PRIM_DISPLAY_P3:return ZIMG_PRIMARIES_ST432_1;
    case MP_CSP_PRIM_APPLE:     // ?
    case MP_CSP_PRIM_ADOBE:
    case MP_CSP_PRIM_PRO_PHOTO:
    case MP_CSP_PRIM_V_GAMUT:
    case MP_CSP_PRIM_S_GAMUT:   // ?
    default:                    return ZIMG_PRIMARIES_BT709;
    }
}

static void destroy_zimg(struct mp_zimg_context *ctx)
{
    talloc_free(ctx->zimg_tmp_alloc);
    ctx->zimg_tmp = ctx->zimg_tmp_alloc = NULL;
    zimg_filter_graph_free(ctx->zimg_graph);
    ctx->zimg_graph = NULL;
    TA_FREEP(&ctx->zimg_src);
    TA_FREEP(&ctx->zimg_dst);
}

static void free_mp_zimg(void *p)
{
    struct mp_zimg_context *ctx = p;

    destroy_zimg(ctx);
}

struct mp_zimg_context *mp_zimg_alloc(void)
{
    struct mp_zimg_context *ctx = talloc_ptrtype(NULL, ctx);
    *ctx = (struct mp_zimg_context) {
        .log = mp_null_log,
    };
    ctx->opts = *(struct zimg_opts *)zimg_conf.defaults;
    talloc_set_destructor(ctx, free_mp_zimg);
    return ctx;
}

void mp_zimg_enable_cmdline_opts(struct mp_zimg_context *ctx,
                                 struct mpv_global *g)
{
    if (ctx->opts_cache)
        return;

    ctx->opts_cache = m_config_cache_alloc(ctx, g, &zimg_conf);
    destroy_zimg(ctx); // force update
    mp_zimg_update_from_cmdline(ctx); // first update
}

static int repack_entrypoint(void *user, unsigned i, unsigned x0, unsigned x1)
{
    struct mp_zimg_repack *r = user;

    // If reading is not aligned, just read slightly more data.
    if (!r->pack)
        x0 &= ~(unsigned)(mp_repack_get_align_x(r->repack) - 1);

    // mp_repack requirements and zimg guarantees.
    assert(!(i & (mp_repack_get_align_y(r->repack) - 1)));
    assert(!(x0 & (mp_repack_get_align_x(r->repack) - 1)));

    unsigned i_src = i & (r->pack ? r->zmask[0] : ZIMG_BUFFER_MAX);
    unsigned i_dst = i & (r->pack ? ZIMG_BUFFER_MAX : r->zmask[0]);

    repack_line(r->repack, x0, i_dst, x0, i_src, x1 - x0);

    return 0;
}

static bool wrap_buffer(struct mp_zimg_repack *r,
                        zimg_image_buffer *buf,
                        struct mp_image *mpi)
{
    *buf = (zimg_image_buffer){ZIMG_API_VERSION};

    bool direct[MP_MAX_PLANES] = {0};

    for (int p = 0; p < mpi->num_planes; p++) {
        // If alignment is good, try to avoid copy.
        direct[p] = !((uintptr_t)mpi->planes[p] % ZIMG_ALIGN) &&
                    !(mpi->stride[p] % ZIMG_ALIGN);
    }

    if (!repack_config_buffers(r->repack, 0, r->pack ? mpi : r->tmp,
                                          0, r->pack ? r->tmp : mpi, direct))
        return false;

    for (int n = 0; n < MP_ARRAY_SIZE(buf->plane); n++) {
        // Note: this is really the only place we have to care about plane
        // permutation (zimg_image_buffer may have a different plane order
        // than the shadow mpi like r->tmp). We never use the zimg indexes
        // in other places.
        int mplane = r->z_planes[n];
        if (mplane < 0)
            continue;

        struct mp_image *tmpi = direct[mplane] ? mpi : r->tmp;
        buf->plane[n].data = tmpi->planes[mplane];
        buf->plane[n].stride = tmpi->stride[mplane];
        buf->plane[n].mask = direct[mplane] ? ZIMG_BUFFER_MAX : r->zmask[mplane];
    }

    return true;
}

// (ctx can be NULL for probing.)
static bool setup_format(zimg_image_format *zfmt, struct mp_zimg_repack *r,
                         bool pack, struct mp_image_params *user_fmt,
                         struct mp_zimg_context *ctx)
{
    r->fmt = *user_fmt;
    r->pack = pack;

    zimg_image_format_default(zfmt, ZIMG_API_VERSION);

    int rp_flags = 0;

    // For e.g. RGB565, go to lowest depth on pack for less weird dithering.
    if (r->pack) {
        rp_flags |= REPACK_CREATE_ROUND_DOWN;
    } else {
        rp_flags |= REPACK_CREATE_EXPAND_8BIT;
    }

    r->repack = mp_repack_create_planar(r->fmt.imgfmt, r->pack, rp_flags);
    if (!r->repack)
        return false;

    int align_x = mp_repack_get_align_x(r->repack);

    r->zimgfmt = r->pack ? mp_repack_get_format_src(r->repack)
                         : mp_repack_get_format_dst(r->repack);

    if (ctx) {
        talloc_steal(r, r->repack);
    } else {
        TA_FREEP(&r->repack);
    }

    struct mp_image_params fmt = r->fmt;
    mp_image_params_guess_csp(&fmt);

    struct mp_regular_imgfmt desc;
    if (!mp_get_regular_imgfmt(&desc, r->zimgfmt))
        return false;

    // Relies on zimg callbacks reading on 64 byte alignment.
    if (!MP_IS_POWER_OF_2(align_x) || align_x > 64 / desc.component_size)
        return false;

    // no weird stuff
    if (desc.num_planes > 4)
        return false;

    for (int n = 0; n < 4; n++)
        r->z_planes[n] = -1;

    for (int n = 0; n < desc.num_planes; n++) {
        if (desc.planes[n].num_components != 1)
            return false;
        int c = desc.planes[n].components[0];
        if (c < 1 || c > 4)
            return false;
        if (c < 4) {
            // Unfortunately, ffmpeg prefers GBR order for planar RGB, while zimg
            // is sane. This makes it necessary to determine and fix the order.
            r->z_planes[c - 1] = n;
        } else {
            r->z_planes[3] = n; // alpha, always plane 4 in zimg

#if HAVE_ZIMG_ALPHA
            zfmt->alpha = fmt.alpha == MP_ALPHA_PREMUL
                ? ZIMG_ALPHA_PREMULTIPLIED : ZIMG_ALPHA_STRAIGHT;
#else
            return false;
#endif
        }
    }

    r->num_planes = desc.num_planes;

    // Note: formats with subsampled chroma may have odd width or height in
    // mpv and FFmpeg. This is because the width/height is actually a cropping
    // rectangle. Reconstruct the image allocation size and set the cropping.
    zfmt->width = r->real_w = MP_ALIGN_UP(fmt.w, 1 << desc.chroma_xs);
    zfmt->height = r->real_h = MP_ALIGN_UP(fmt.h, 1 << desc.chroma_ys);
    if (!r->pack && ctx) {
        // Relies on ctx->zimg_dst being initialized first.
        struct mp_zimg_repack *dst = ctx->zimg_dst;
        zfmt->active_region.width = dst->real_w * (double)fmt.w / dst->fmt.w;
        zfmt->active_region.height = dst->real_h * (double)fmt.h / dst->fmt.h;
    }

    zfmt->subsample_w = desc.chroma_xs;
    zfmt->subsample_h = desc.chroma_ys;

    zfmt->color_family = ZIMG_COLOR_YUV;
    if (desc.num_planes <= 2) {
        zfmt->color_family = ZIMG_COLOR_GREY;
    } else if (fmt.color.space == MP_CSP_RGB || fmt.color.space == MP_CSP_XYZ) {
        zfmt->color_family = ZIMG_COLOR_RGB;
    }

    if (desc.component_type == MP_COMPONENT_TYPE_UINT &&
        desc.component_size == 1)
    {
        zfmt->pixel_type = ZIMG_PIXEL_BYTE;
    } else if (desc.component_type == MP_COMPONENT_TYPE_UINT &&
               desc.component_size == 2)
    {
        zfmt->pixel_type = ZIMG_PIXEL_WORD;
    } else if (desc.component_type == MP_COMPONENT_TYPE_FLOAT &&
               desc.component_size == 2)
    {
        zfmt->pixel_type = ZIMG_PIXEL_HALF;
    } else if (desc.component_type == MP_COMPONENT_TYPE_FLOAT &&
               desc.component_size == 4)
    {
        zfmt->pixel_type = ZIMG_PIXEL_FLOAT;
    } else {
        return false;
    }

    // (Formats like P010 are basically reported as P016.)
    zfmt->depth = desc.component_size * 8 + MPMIN(0, desc.component_pad);

    zfmt->pixel_range = fmt.color.levels == MP_CSP_LEVELS_PC ?
                        ZIMG_RANGE_FULL : ZIMG_RANGE_LIMITED;

    zfmt->matrix_coefficients = mp_to_z_matrix(fmt.color.space);
    zfmt->transfer_characteristics = mp_to_z_trc(fmt.color.gamma);
    zfmt->color_primaries = mp_to_z_prim(fmt.color.primaries);
    zfmt->chroma_location = mp_to_z_chroma(fmt.chroma_location);

    if (ctx && ctx->opts.fast) {
        // mpv's default for RGB output slows down zimg significantly.
        if (zfmt->transfer_characteristics == ZIMG_TRANSFER_IEC_61966_2_1 &&
            zfmt->color_family == ZIMG_COLOR_RGB)
            zfmt->transfer_characteristics = ZIMG_TRANSFER_BT709;
    }

    // mpv treats _some_ gray formats as RGB; zimg doesn't like this.
    if (zfmt->color_family == ZIMG_COLOR_GREY &&
        zfmt->matrix_coefficients == ZIMG_MATRIX_RGB)
        zfmt->matrix_coefficients = ZIMG_MATRIX_BT470_BG;

    return true;
}

static bool allocate_buffer(struct mp_zimg_context *ctx,
                            struct mp_zimg_repack *r)
{
    unsigned lines = 0;
    int err;
    if (r->pack) {
        err = zimg_filter_graph_get_output_buffering(ctx->zimg_graph, &lines);
    } else {
        err = zimg_filter_graph_get_input_buffering(ctx->zimg_graph, &lines);
    }

    if (err)
        return false;

    r->zmask[0] = zimg_select_buffer_mask(lines);

    // Either ZIMG_BUFFER_MAX, or a power-of-2 slice buffer.
    assert(r->zmask[0] == ZIMG_BUFFER_MAX || MP_IS_POWER_OF_2(r->zmask[0] + 1));

    int h = r->zmask[0] == ZIMG_BUFFER_MAX ? r->fmt.h : r->zmask[0] + 1;
    if (h >= r->fmt.h) {
        h = r->fmt.h;
        r->zmask[0] = ZIMG_BUFFER_MAX;
    }

    r->tmp = mp_image_alloc(r->zimgfmt, r->fmt.w, h);
    talloc_steal(r, r->tmp);

    if (!r->tmp)
        return false;

    // Note: although zimg doesn't require that the chroma plane's zmask is
    //       divided by the full size zmask, the repack callback requires it,
    //       since mp_repack can handle only proper slices.
    for (int n = 1; n < r->tmp->fmt.num_planes; n++) {
        r->zmask[n] = r->zmask[0];
        if (r->zmask[0] != ZIMG_BUFFER_MAX)
            r->zmask[n] = r->zmask[n] >> r->tmp->fmt.ys[n];
    }

    return true;
}

bool mp_zimg_config(struct mp_zimg_context *ctx)
{
    struct zimg_opts *opts = &ctx->opts;

    destroy_zimg(ctx);

    if (ctx->opts_cache)
        mp_zimg_update_from_cmdline(ctx);

    ctx->zimg_src = talloc_zero(NULL, struct mp_zimg_repack);
    ctx->zimg_dst = talloc_zero(NULL, struct mp_zimg_repack);

    zimg_image_format src_fmt, dst_fmt;

    // Note: do zimg_dst first, because zimg_src uses fields from zimg_dst.
    if (!setup_format(&dst_fmt, ctx->zimg_dst, true, &ctx->dst, ctx) ||
        !setup_format(&src_fmt, ctx->zimg_src, false, &ctx->src, ctx))
        goto fail;

    zimg_graph_builder_params params;
    zimg_graph_builder_params_default(&params, ZIMG_API_VERSION);

    params.resample_filter = opts->scaler;
    params.filter_param_a = opts->scaler_params[0];
    params.filter_param_b = opts->scaler_params[1];

    params.resample_filter_uv = opts->scaler_chroma;
    params.filter_param_a_uv = opts->scaler_chroma_params[0];
    params.filter_param_b_uv = opts->scaler_chroma_params[1];

    params.dither_type = opts->dither;

    params.cpu_type = ZIMG_CPU_AUTO_64B;

    if (opts->fast)
        params.allow_approximate_gamma = 1;

    if (ctx->src.color.sig_peak > 0)
        params.nominal_peak_luminance = ctx->src.color.sig_peak;

    ctx->zimg_graph = zimg_filter_graph_build(&src_fmt, &dst_fmt, &params);
    if (!ctx->zimg_graph) {
        char err[128] = {0};
        zimg_get_last_error(err, sizeof(err) - 1);
        MP_ERR(ctx, "zimg_filter_graph_build: %s \n", err);
        goto fail;
    }

    size_t tmp_size;
    if (!zimg_filter_graph_get_tmp_size(ctx->zimg_graph, &tmp_size)) {
        tmp_size = MP_ALIGN_UP(tmp_size, ZIMG_ALIGN) + ZIMG_ALIGN;
        ctx->zimg_tmp_alloc = ta_alloc_size(NULL, tmp_size);
        if (ctx->zimg_tmp_alloc) {
            ctx->zimg_tmp =
                (void *)MP_ALIGN_UP((uintptr_t)ctx->zimg_tmp_alloc, ZIMG_ALIGN);
        }
    }

    if (!ctx->zimg_tmp_alloc)
        goto fail;

    if (!allocate_buffer(ctx, ctx->zimg_src) ||
        !allocate_buffer(ctx, ctx->zimg_dst))
        goto fail;

    return true;

fail:
    destroy_zimg(ctx);
    return false;
}

bool mp_zimg_config_image_params(struct mp_zimg_context *ctx)
{
    if (ctx->zimg_src && mp_image_params_equal(&ctx->src, &ctx->zimg_src->fmt) &&
        ctx->zimg_dst && mp_image_params_equal(&ctx->dst, &ctx->zimg_dst->fmt) &&
        (!ctx->opts_cache || !m_config_cache_update(ctx->opts_cache)) &&
        ctx->zimg_graph)
        return true;
    return mp_zimg_config(ctx);
}

bool mp_zimg_convert(struct mp_zimg_context *ctx, struct mp_image *dst,
                     struct mp_image *src)
{
    ctx->src = src->params;
    ctx->dst = dst->params;

    if (!mp_zimg_config_image_params(ctx)) {
        MP_ERR(ctx, "zimg initialization failed.\n");
        return false;
    }

    assert(ctx->zimg_graph);

    zimg_image_buffer zsrc, zdst;
    if (!wrap_buffer(ctx->zimg_src, &zsrc, src) ||
        !wrap_buffer(ctx->zimg_dst, &zdst, dst))
    {
        MP_ERR(ctx, "zimg repacker initialization failed.\n");
        return false;
    }

    // An annoyance.
    zimg_image_buffer_const zsrc_c = {ZIMG_API_VERSION};
    for (int n = 0; n < MP_ARRAY_SIZE(zsrc_c.plane); n++) {
        zsrc_c.plane[n].data = zsrc.plane[n].data;
        zsrc_c.plane[n].stride = zsrc.plane[n].stride;
        zsrc_c.plane[n].mask = zsrc.plane[n].mask;
    }

    // (The API promises to succeed if no user callbacks fail, so no need
    // to check the return value.)
    zimg_filter_graph_process(ctx->zimg_graph, &zsrc_c, &zdst,
                              ctx->zimg_tmp,
                              repack_entrypoint, ctx->zimg_src,
                              repack_entrypoint, ctx->zimg_dst);

    return true;
}

static bool supports_format(int imgfmt, bool out)
{
    struct mp_image_params fmt = {.imgfmt = imgfmt};
    struct mp_zimg_repack t;
    zimg_image_format zfmt;
    return setup_format(&zfmt, &t, out, &fmt, NULL);
}

bool mp_zimg_supports_in_format(int imgfmt)
{
    return supports_format(imgfmt, false);
}

bool mp_zimg_supports_out_format(int imgfmt)
{
    return supports_format(imgfmt, true);
}
