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

#include <stddef.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>
#include <inttypes.h>

#include "common/common.h"
#include "draw_bmp.h"
#include "img_convert.h"
#include "video/mp_image.h"
#include "video/repack.h"
#include "video/sws_utils.h"
#include "video/img_format.h"
#include "video/csputils.h"

const bool mp_draw_sub_formats[SUBBITMAP_COUNT] = {
    [SUBBITMAP_LIBASS] = true,
    [SUBBITMAP_RGBA] = true,
};

struct part {
    int change_id;
    // Sub-bitmaps scaled to final sizes.
    int num_imgs;
    struct mp_image **imgs;
};

// Must be powers of 2.
#define TILE_W 256
#define TILE_H 64

struct extent {
    uint16_t x0, x1;
};

struct tile {
    // X range of non-transparent pixels for each line in this tile.
    struct extent exts[TILE_H];
    bool any_pixels;
};

struct mp_draw_sub_cache
{
    // Possibly cached parts. Also implies what's in the video_overlay.
    struct part parts[MAX_OSD_PARTS];

    struct mp_image_params params;  // target image params

    int align_x, align_y;           // alignment for all video pixels

    struct mp_image *rgba_overlay;  // all OSD in RGBA
    struct mp_image *video_overlay; // rgba_overlay converted to video colorspace

    unsigned t_w, t_h;              // size in tiles (bottom/right pixels cut off)
    struct tile *tiles;             // tiles[y_tile * t_w + x_tile]

    struct mp_sws_context *rgba_to_overlay; // scaler for rgba -> video csp.

    struct mp_repack *overlay_to_f32; // convert video_overlay to float
    struct mp_image *overlay_tmp;   // slice in float32

    struct mp_repack *video_to_f32; // convert video to float
    struct mp_repack *video_from_f32; // convert float back to video
    struct mp_image *video_tmp;     // slice in float32

    // TODO: add a gray scaler to do the cursed alpha->chroma down sampling
};

static void blend_slice(struct mp_draw_sub_cache *cache, int rgb_y)
{
    struct mp_image *ov = cache->overlay_tmp;
    struct mp_image *vid = cache->video_tmp;

    ptrdiff_t astride = ov->stride[ov->num_planes - 1] / 4; // don't do this at home

    uint32_t *rgba = mp_image_pixel_ptr(cache->rgba_overlay, 0, 0, rgb_y);

    for (int p = 0; p < vid->num_planes; p++) {
        int xs = vid->fmt.xs[p];
        int ys = vid->fmt.ys[p];
        int h = (1 << vid->fmt.chroma_ys) - (1 << ys) + 1;
        int cw = mp_chroma_div_up(vid->w, xs);
        for (int y = 0; y < h; y++) {
            float *ov_ptr = mp_image_pixel_ptr(ov, p, 0, y);
            float *a_ptr = mp_image_pixel_ptr(ov, ov->num_planes - 1, 0, y);
            float *vid_ptr = mp_image_pixel_ptr(vid, p, 0, y);

            for (int x = 0; x < cw; x++) {
                // hurrrrrrrrrrrr
                float alpha = 0;
                for (int ax = 0; ax < (1 << xs); ax++) {
                    for (int ay = 0; ay < (1 << ys); ay++) {
                        alpha += a_ptr[(x << xs) + ax + astride * ay];
                    }
                }
                alpha /= (1 << xs) * (1 << ys);

                vid_ptr[x] = ov_ptr[x] + (1 - alpha) * vid_ptr[x];
            }
        }
    }
}

static void blend_overlay_with_video(struct mp_draw_sub_cache *cache,
                                     struct mp_image *dst)
{
    if (!repack_config_buffers(cache->video_to_f32,
                               0, cache->video_tmp, 0, dst, NULL))
        printf("fuck! 1\n");
    if (!repack_config_buffers(cache->video_from_f32,
                               0, dst, 0, cache->video_tmp, NULL))
        printf("fuck! 2\n");

    // TODO: use extents to reduce blending

    for (int y = 0; y < dst->h; y += cache->align_y) {
        repack_line(cache->overlay_to_f32, 0, 0, 0, y, dst->w);
        repack_line(cache->video_to_f32, 0, 0, 0, y, dst->w);

        blend_slice(cache, y);

        repack_line(cache->video_from_f32, 0, y, 0, 0, dst->w);
    }
}

static void convert_to_video_overlay(struct mp_draw_sub_cache *cache)
{
    if (!cache->video_overlay)
        return;

    // TODO: use tile shit to avoid converting some parts of the overlay
    if (mp_sws_scale(cache->rgba_to_overlay, cache->video_overlay,
                     cache->rgba_overlay) < 0)
        printf("eh what\n");
}

static void draw_ass_rgba(unsigned char *src, int src_w, int src_h,
                          int src_stride, unsigned char *dst, size_t dst_stride,
                          int dst_x, int dst_y, uint32_t color)
{
    const unsigned int r = (color >> 24) & 0xff;
    const unsigned int g = (color >> 16) & 0xff;
    const unsigned int b = (color >>  8) & 0xff;
    const unsigned int a = 0xff - (color & 0xff);

    dst += dst_y * dst_stride + dst_x * 4;

    for (int y = 0; y < src_h; y++, dst += dst_stride, src += src_stride) {
        uint32_t *dstrow = (uint32_t *) dst;
        for (int x = 0; x < src_w; x++) {
            const unsigned int v = src[x];
            int rr = (r * a * v);
            int gg = (g * a * v);
            int bb = (b * a * v);
            int aa =      a * v;
            uint32_t dstpix = dstrow[x];
            unsigned int dstb =  dstpix        & 0xFF;
            unsigned int dstg = (dstpix >>  8) & 0xFF;
            unsigned int dstr = (dstpix >> 16) & 0xFF;
            unsigned int dsta = (dstpix >> 24) & 0xFF;
            dstb = (bb       + dstb * (255 * 255 - aa)) / (255 * 255);
            dstg = (gg       + dstg * (255 * 255 - aa)) / (255 * 255);
            dstr = (rr       + dstr * (255 * 255 - aa)) / (255 * 255);
            dsta = (aa * 255 + dsta * (255 * 255 - aa)) / (255 * 255);
            dstrow[x] = dstb | (dstg << 8) | (dstr << 16) | (dsta << 24);
        }
    }
}

static void render_ass(struct mp_draw_sub_cache *cache, struct sub_bitmaps *sb)
{
    assert(sb->format == SUBBITMAP_LIBASS);

    // TODO: write extents

    for (int i = 0; i < sb->num_parts; i++) {
        struct sub_bitmap *s = &sb->parts[i];

        draw_ass_rgba(s->bitmap, s->w, s->h, s->stride,
                      cache->rgba_overlay->planes[0],
                      cache->rgba_overlay->stride[0],
                      s->x, s->y, s->libass.color);
    }
}

static void render_sb(struct mp_draw_sub_cache *cache, struct sub_bitmaps *sb)
{
    struct part *part = &cache->parts[sb->render_index];

    if (sb->format == SUBBITMAP_LIBASS) {
        render_ass(cache, sb);
    } else {
        printf("TODO: rgba sub-bitmaps unimplemented\n");
    }

    part->change_id = sb->change_id;
}

static void clear_rgba_overlay(struct mp_draw_sub_cache *cache)
{
    // TODO: use extents crap to clear only uncleared areas?

    assert(cache->rgba_overlay->imgfmt == IMGFMT_BGR32);
    for (int y = 0; y < cache->rgba_overlay->h; y++) {
        memset(mp_image_pixel_ptr(cache->rgba_overlay, 0, 0, y), 0,
               cache->rgba_overlay->w * 4);
    }
}

static bool reinit(struct mp_draw_sub_cache *cache, struct mp_image_params *params)
{
    talloc_free_children(cache);
    *cache = (struct mp_draw_sub_cache){.params = *params};

    int rflags = REPACK_CREATE_EXPAND_8BIT | REPACK_CREATE_PLANAR_F32;

    cache->video_to_f32 = mp_repack_create_planar(params->imgfmt, false, rflags);
    talloc_steal(cache, cache->video_to_f32);
    if (!cache->video_to_f32)
        return false;

    int vid_f32_fmt = mp_repack_get_format_dst(cache->video_to_f32);

    cache->video_from_f32 = mp_repack_create_planar(params->imgfmt, true, rflags);
    talloc_steal(cache, cache->video_from_f32);
    if (!cache->video_from_f32)
        return false;

    assert(mp_repack_get_format_dst(cache->video_to_f32) ==
           mp_repack_get_format_src(cache->video_from_f32));

    // Find a reasonable intermediate format for video_overlay. Requirements:
    //  - same subsampling
    //  - has alpha
    //  - uses video colorspace
    //  - REPACK_CREATE_PLANAR_F32 support
    //  - converted float values for video and overlay use same value ranges
    //  - probably not using float (vaguely wastes memory)
    struct mp_regular_imgfmt vfdesc = {0};
    mp_get_regular_imgfmt(&vfdesc, mp_repack_get_format_dst(cache->video_to_f32));
    assert(vfdesc.component_type == MP_COMPONENT_TYPE_FLOAT);
    struct mp_regular_imgfmt odesc = vfdesc;

    int overlay_fmt = 0;
    if (params->color.space == MP_CSP_RGB) {
        // No point in doing anything fancy.
        // TODO: check for more corner cases (different gamma/primaries?)
        overlay_fmt = IMGFMT_BGR32;
    } else {
        // Try to use video component type to reduce problems with value ranges.
        // The low bit depth ones typically don't have alpha variants either.
        struct mp_regular_imgfmt vdesc = {0};
        int vdepth = 0;
        if (mp_get_regular_imgfmt(&vdesc, params->imgfmt))
            vdepth = vdesc.component_size * 8 + MPMIN(vdesc.component_pad, 0);
        if (vdepth >= 8) {
            odesc.component_type = vdesc.component_type;
            odesc.component_size = vdesc.component_size;
            odesc.component_pad = MPMIN(vdesc.component_pad, 0);
        } else {
            odesc.component_type = MP_COMPONENT_TYPE_UINT;
            odesc.component_size = 1;
            odesc.component_pad = 0;
        }

        // Ensure there's alpha.
        if (odesc.planes[odesc.num_planes - 1].components[0] != 4) {
            if (odesc.num_planes >= 4)
                return false; // wat
            odesc.planes[odesc.num_planes++] =
                (struct mp_regular_imgfmt_plane){1, {4}};
        }

        overlay_fmt = mp_find_regular_imgfmt(&odesc);
    }
    if (!overlay_fmt)
        return false;

    cache->overlay_to_f32 = mp_repack_create_planar(overlay_fmt, false, rflags);
    talloc_steal(cache, cache->overlay_to_f32);
    if (!cache->overlay_to_f32)
        return false;

    int render_fmt = mp_repack_get_format_dst(cache->overlay_to_f32);

    struct mp_regular_imgfmt ofdesc = {0};
    mp_get_regular_imgfmt(&ofdesc, render_fmt);

    if (ofdesc.planes[ofdesc.num_planes - 1].components[0] != 4)
        return false;

    // The formats must be the same, minus possible lack of alpha in vfdesc.
    if (ofdesc.num_planes != vfdesc.num_planes &&
        ofdesc.num_planes - 1 != vfdesc.num_planes)
        return false;
    for (int n = 0; n < vfdesc.num_planes; n++) {
        if (vfdesc.planes[n].components[0] != ofdesc.planes[n].components[0])
            return false;
    }

    cache->align_x = MPMAX(mp_repack_get_align_x(cache->overlay_to_f32),
                           mp_repack_get_align_x(cache->video_to_f32));
    cache->align_y = MPMAX(mp_repack_get_align_y(cache->overlay_to_f32),
                           mp_repack_get_align_y(cache->video_to_f32));

    assert(cache->align_x < TILE_W && cache->align_y < TILE_H);

    int w = MP_ALIGN_UP(params->w, cache->align_x);
    int h = MP_ALIGN_UP(params->h, cache->align_y);

    cache->rgba_overlay = talloc_steal(cache, mp_image_alloc(IMGFMT_BGR32, w, h));
    cache->overlay_tmp = talloc_steal(cache,
                            mp_image_alloc(render_fmt, w /*TILE_W*/, cache->align_y));
    cache->video_tmp = talloc_steal(cache,
                            mp_image_alloc(vid_f32_fmt, w /*TILE_W*/, cache->align_y));
    if (!cache->rgba_overlay || !cache->overlay_tmp || !cache->video_tmp)
        return false;

    mp_image_params_guess_csp(&cache->rgba_overlay->params);
    cache->rgba_overlay->params.alpha = MP_ALPHA_PREMUL;

    cache->overlay_tmp->params.color = params->color;
    cache->video_tmp->params.color = params->color;

    if (cache->rgba_overlay->imgfmt == overlay_fmt) {
        if (!repack_config_buffers(cache->overlay_to_f32, 0, cache->overlay_tmp,
                                   0, cache->rgba_overlay, NULL))
        return false;
    } else {
        cache->video_overlay = talloc_steal(cache,
                                            mp_image_alloc(overlay_fmt, w, h));
        if (!cache->video_overlay)
            return false;

        cache->video_overlay->params.color = params->color;
        cache->video_overlay->params.chroma_location = params->chroma_location;
        cache->video_overlay->params.alpha = MP_ALPHA_PREMUL;

        cache->rgba_to_overlay = mp_sws_alloc(cache);
        cache->rgba_to_overlay->force_scaler = MP_SWS_ZIMG;
        if (!mp_sws_supports_formats(cache->rgba_to_overlay,
                    cache->video_overlay->imgfmt, cache->rgba_overlay->imgfmt))
            return false;

        if (!repack_config_buffers(cache->overlay_to_f32, 0, cache->overlay_tmp,
                                   0, cache->video_overlay, NULL))
            return false;
    }

    cache->t_w = MP_ALIGN_UP(cache->rgba_overlay->w, TILE_W) / TILE_W;
    cache->t_h = MP_ALIGN_UP(cache->rgba_overlay->h, TILE_H) / TILE_H;

    cache->tiles = talloc_zero_array(cache, struct tile, cache->t_w * cache->t_h);

    printf("ov_vid: %s\n", vo_format_name(cache->video_overlay ? cache->video_overlay->imgfmt : 0));
    printf("ov_f32: %s\n", vo_format_name(cache->overlay_tmp->imgfmt));
    printf("vid_f32: %s\n", vo_format_name(cache->video_tmp->imgfmt));
    printf("align=%d:%d\n", cache->align_x, cache->align_y);

    // TODO: init extents

    return true;
}

// cache: if not NULL, the function will set *cache to a talloc-allocated cache
//        containing scaled versions of sbs contents - free the cache with
//        talloc_free()
void mp_draw_sub_bitmaps(struct mp_draw_sub_cache **p_cache, struct mp_image *dst,
                         struct sub_bitmap_list *sbs_list)
{
    struct mp_draw_sub_cache *cache = p_cache ? *p_cache : NULL;
    if (!cache)
        cache = talloc_zero(NULL, struct mp_draw_sub_cache);

    if (!mp_image_params_equal(&cache->params, &dst->params) || !cache->video_tmp)
    {
        if (!reinit(cache, &dst->params)) {
            talloc_free_children(cache);
            *cache = (struct mp_draw_sub_cache){0};
            printf("failed\n");
            goto done;
        }
    }

    struct sub_bitmaps *parts[MAX_OSD_PARTS] = {0};
    bool dirty = false;

    for (int n = 0; n < sbs_list->num_items; n++)
        parts[sbs_list->items[n]->render_index] = sbs_list->items[n];

    for (int n = 0; n < MAX_OSD_PARTS; n++) {
        struct sub_bitmaps *sb = parts[n];
        int change_id = sb ? sb->change_id : 0;
        dirty |= change_id != cache->parts[n].change_id;
    }

    if (dirty) {
        clear_rgba_overlay(cache);

        for (int n = 0; n < MAX_OSD_PARTS; n++) {
            struct sub_bitmaps *sb = parts[n];
            if (sb) {
                render_sb(cache, sb);
            } else {
                cache->parts[n].change_id = 0;
            }
        }

        convert_to_video_overlay(cache);
    }

    blend_overlay_with_video(cache, dst);

done:
    if (p_cache) {
        *p_cache = cache;
    } else {
        talloc_free(cache);
    }
}

// vim: ts=4 sw=4 et tw=80
