/*
    3D Lena heightmap demo.
    Decodes the compressed Lena image from externals/lena.c, lifts luminance
    into a heightfield, and renders colored quads via the public b3d API.
    Supports headless PNG snapshots with --snapshot=PATH or B3D_SNAPSHOT.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "SDL.h"
#include "b3d.h"
#include "pngwrite.h"

/* ------------------------- PNG helpers -------------------------- */
static const char *get_snapshot_path(int argc, char **argv)
{
    const char *env = getenv("B3D_SNAPSHOT");
    if (env && env[0]) return env;
    for (int i = 1; i < argc; ++i) {
        const char *flag = "--snapshot=";
        size_t len = strlen(flag);
        if (strncmp(argv[i], flag, len) == 0) {
            return argv[i] + (int)len;
        }
    }
    return NULL;
}

static void write_png(const char *path, const uint32_t *rgba, int width, int height)
{
    FILE *file = fopen(path, "wb");
    if (!file) return;
    uint8_t *out = malloc((size_t)width * (size_t)height * 4);
    if (!out) {
        fclose(file);
        return;
    }
    size_t idx = 0;
    for (int i = 0; i < width * height; ++i) {
        uint32_t p = rgba[i];
        out[idx++] = (uint8_t)((p >> 16) & 0xff);
        out[idx++] = (uint8_t)((p >> 8) & 0xff);
        out[idx++] = (uint8_t)(p & 0xff);
        out[idx++] = 0xff;
    }
    png_write(file, (unsigned)width, (unsigned)height, out, true);
    free(out);
    fclose(file);
}

/* ----------------------- Lena decoder --------------------------- */

#define ACTX_SIGN 3
#define ACTX_VDATA 4
#define ACTX_LEN 5
#define ACTX_LEVEL 25
#define ACTX_IPRED 73
#define ACTX_UE_LEN 10
#define ACTX_COUNT2 166
#define ACTX_EOB2 61

static char *inp =
    "	{k/;	y{ q ; }	c {	@;	={ 	S}	c}	W;;	{4}	k "
    "|; w{	+9;{;	8; 9{	S;	/}	y{ K}	{;}	l{	{ ~{ ;	V}"
    "k}g< t{	E	v;M{ B}y}	<{7;/;	Y} t}kp; Y} $Ha{e} "
    "w};} R} /{>}a	;} ;	`	$W-}	D}B; e;f;*;	~;A;s "
    "O{	o;>{1; m{ `} R}]{ T} v}={ I} ; }a?&; A}$;W;R{u} `; j}W;"
    "s{e}	A;[	R;	X  P; 4 ,F;({<8{#;%}@J{)}	}o^*{u/{"
    "'}]{	*}	}	;{ r}	f	/;}e} }w{ ${{;,; @ d	$}];"
    ">(}	I{ d}	&;	U}	{	y;Y}	{ P{	R} T}_{ }R } l	{ T}"
    "';	|; ${=}	H} (}}8{cp{ s} #}+}	3}kF}<H	 .{ }G}"
    "x;	r	D c{; W; {	b;6; k{}B;*};	]} ~	{ ;;} !}}	x}"
    "v}n;^;	6V}Y{ h; ~	%*}! H; G{ r{ f;Y{ i}z} N  %}.{;	( "
    "	v} _}	h; 7;<}	^;Z;0; ;	<;<; M; N{	}	_{O} !{f{]{"
    "M{;A{}	0;S}${	@;x}y}@	L;1	t{ 3{c{s{_{	`{	D{ ]}"
    "!;	${	_J;v+ }	3{B; ]{	}	E6	.x{?+; {x; }v{$};6}T; "
    "O; ; (}X7}	j; @} :}#	c{ !{ }x	KXt} >; ?{ c; ;	W;	; l;} "
    "h}p}	i{ %	}P}	/{	*}	%L; ;	!{ S{ n} "
    "x;  { 1	J;v{	U}({	@ X{ k} H;4;e J	6;;v; G{{]	&{"
    "A d{ lM{;K;;	4-{}} p h{;	{	rW;	v{;	f}	}1{^&{9{"
    "{ ;~;n;q{	9 R	6{	{ u;a;	;	U;	;Y}	+}}2sk; 8	{	J"
    "K;'i;	;$;	W{	P!{{{P	} [;	(;Q; Un;+}g{C;{"
    "{	; <{	vS} b;6`} ?{+	%;	}n;q{ r}k; ;{c{ S} 2}"
    "~{	4;RW v} R;	kI}|; d; [ O}5; ;;}Z d	{ {&;h	o{ "
    "V	v ;	_{{/}  F{f{r{4{{?{ 4;S}	:;];E}	;	&} #e !{"
    ">{H; {O{ 0;} H;	p; w}>{1}{	-} 4;"
    "S}}	u L{ y} %;2  |{(}	/;,{ )}Y;g}	G}v;T}	};}i {{"
    "};[{ E{q} g;T{ ={}R;	k{ j;_;h}gPc;({	F;6}	}} 3	,}<; "
    "0	 P;{'t}u};		}U}s{8{ E} >{}E	{G{H :{  Yo"
    "g}	}F  D{ R{	 -;M?;= q}_ U	{ ;	 I	{ |{{}	 	1{"
    ",}{ x{{ U{ s;J}}	6{>7;,{ D{	{{ ;]}	;M; &}{ V}	"
    "n{&	T~;({	}[;	r{#	u{X 9;L; Uf})}   {T}		p{	N;	"
    ">{	>	}}D} m{1{	{}X; o}	w}$}	^v} K  f	,}	^3; "
    "{ @{_} _{	o;	4}	h}H;#.{	{}	;	<{ {G{ $;{ "
    "z {a{{D;	?|}{{ ;	`} }	Q}j;4} 	3{Q}	{	* ;}r{"
    "a}	} R{p @;  N{ {f; A;8}L	$}{ }}J{ }	k{r} { [; "
    "-;p{	I{ {	&}J;	T}	?{Z{>;	5>; ];  wz ^}	u;);	H}	; "
    "L	&;	V	E{1{g;C} V} ~;U; ^{	J; { /}	{;(}y} aK /}	.}"
    ";K;N{w{ `{	}T{l`; #;N{lX;	?; +}{ 	w{	;	q;	z;_;"
    "y} 8} 	&{X}	V{ WG}	,; [}U{	v{	Q;	w{	[	Y}N	Yu i{ "
    "{!A{}{ b0;	X~} ;-; 8{	E }	;F{	y{}{	";

#define IMG_SIZE_MAX_LOG2 20
#define DCT_BITS 10
#define DCT_SIZE_LOG2_MAX 5
#define DCT_SIZE_MAX 32
#define DCT_SIZE_MAX4 128
#define DCT_SIZE_MAX_SQ2 2048
#define FREQ_MAX 63
#define SYM_COUNT 1968

static int img_data[3][1 << IMG_SIZE_MAX_LOG2];
static int a_ctx[ACTX_COUNT2];
static int a_low, a_range = 1, stride, y_scale, c_scale;
static int dct_coef[DCT_SIZE_MAX4];

static int get_bit(int c)
{
    int v, *p = a_ctx + c * 2, b = *p + 1, s = b + p[1] + 1;
    if (a_range < SYM_COUNT) {
        a_range *= SYM_COUNT;
        a_low *= SYM_COUNT;
        if ((v = *inp)) {
            a_low += (v - 1 - (v > 10) - (v > 13) - (v > 34) - (v > 92)) << 4;
            v = *++inp;
            inp++;
            a_low += v < 33 ? (v ^ 8) * 2 % 5
                            : (v ^ 6) % 3 * 4 + (*inp++ ^ 8) * 2 % 5 + 4;
        }
    }
    v = a_range * b / s;
    if ((b = (a_low >= v))) {
        a_low -= v;
        a_range -= v;
    } else {
        a_range = v;
    }
    p[b]++;
    if (s > FREQ_MAX) {
        *p /= 2;
        p[1] /= 2;
    }
    return b;
}

static int get_ue(int c)
{
    int i = 0, v = 1;
    while (!get_bit(c + i))
        i++;
    while (i--)
        v += v + get_bit(ACTX_VDATA);
    return v - 1;
}

static void idct(int *dst,
                 int dst_stride,
                 int *src,
                 int src_stride,
                 int stride2,
                 int n,
                 int rshift)
{
    for (int l = 0; l < n; l++)
        for (int i = 0; i < n; i++) {
            int sum = 1 << (rshift - 1);
            for (int j = 0; j < n; j++)
                sum += src[j * src_stride + l * stride2] *
                       dct_coef[(2 * i + 1) * j * DCT_SIZE_MAX / n %
                                DCT_SIZE_MAX4];
            dst[i * dst_stride + l * stride2] = sum >> rshift;
        }
}

static int buf1[DCT_SIZE_MAX_SQ2];
static void decode_rec(int x, int y, int w_log2)
{
    int b;
    int w = 1 << w_log2, n = w * w;

    if ((w_log2 > DCT_SIZE_LOG2_MAX) || (w_log2 > 2 && get_bit(w_log2 - 3))) {
        w /= 2;
        for (int i = 0; i < 4; i++)
            decode_rec(x + i % 2 * w, y + i / 2 * w, w_log2 - 1);
        return;
    }

    int pred_idx = get_ue(ACTX_IPRED);
    for (int c_idx = 0; c_idx < 3; c_idx++) {
        int *out = img_data[c_idx] + y * stride + x;
        int c_idx1 = c_idx > 0;

        memset(buf1, 0, n * sizeof(int));
        for (int i = 0; i < n; i++) {
            if (get_bit(ACTX_EOB2 + w_log2 * 2 + c_idx1))
                break;
            i += get_ue(ACTX_LEN + c_idx1 * ACTX_UE_LEN);
            b = 1 - 2 * get_bit(ACTX_SIGN);
            buf1[i] =
                b *
                (get_ue(ACTX_LEVEL + (c_idx1 + (i < n / 8) * 2) * ACTX_UE_LEN) +
                 1) *
                (c_idx ? c_scale : y_scale);
        }

        if (!pred_idx) {
            int dc = 0;
            for (int i = 0; i < w; i++) {
                dc += y ? out[-stride + i] : 0;
                dc += x ? out[i * stride - 1] : 0;
            }
            *buf1 += x && y ? dc / 2 : dc;
        }

        idct(buf1 + n, 1, buf1, 1, w, w, DCT_BITS);
        idct(out, stride, buf1 + n, w, 1, w, DCT_BITS + w_log2);

        if (!pred_idx)
            continue;

        int swap = pred_idx < 17, frac;
        int delta = swap ? 9 - pred_idx : pred_idx - 25;
        for (int i = 0; i < w; i++)
            for (int j = 0; j < w; j++) {
                for (int k = 0; k < 2; k++) {
                    int x1 = i * delta + delta;
                    frac = x1 & 7;
                    x1 = (x1 >> 3) + j + k;
                    if ((b = (x1 < 0)))
                        x1 = (x1 * 8 + delta / 2) / delta - 2;
                    x1 = x1 < w ? x1 : w - 1;
                    buf1[k] =
                        b ^ swap ? out[x1 * stride - 1] : out[-stride + x1];
                }
                out[swap ? j * stride + i : i * stride + j] +=
                    (*buf1 * (8 - frac) + buf1[1] * frac + 4) >> 3;
            }
    }
}

static uint32_t *decode_lena(int *out_w, int *out_h)
{
    int a = 0;
    int b = 74509276;
    for (int i = 0; i < 128; i++) {
        dct_coef[i + 96 & 127] = ((a >> 19) + 1) >> 1;
        int c = b;
        b = (2144896910LL * b >> 30) - a;
        a = c;
    }

    *dct_coef = 1024;
    int w_log2 = get_ue(ACTX_LEN);
    stride = 1 << w_log2;
    int h = stride - get_ue(ACTX_LEN);
    y_scale = get_ue(ACTX_LEN);
    c_scale = get_ue(ACTX_LEN);

    decode_rec(0, 0, w_log2);

    *out_w = stride;
    *out_h = h;
    size_t count = (size_t)stride * (size_t)h;
    uint32_t *rgba = malloc(count * sizeof(uint32_t));
    if (!rgba)
        return NULL;

    for (int i = 0; i < h * stride; i++) {
        int y = img_data[0][i], cg = img_data[1][i], co = img_data[2][i];
        int t = y - cg;
        int r = t + co;
        int g = y + cg;
        int bcol = t - co;
        if (r < 0) r = 0; else if (r > 255) r = 255;
        if (g < 0) g = 0; else if (g > 255) g = 255;
        if (bcol < 0) bcol = 0; else if (bcol > 255) bcol = 255;
        rgba[i] = (uint32_t)(r << 16 | g << 8 | bcol);
    }
    return rgba;
}

/* ----------------------- Rendering ------------------------------ */

static void render_frame(uint32_t *pixels, b3d_depth_t *depth,
                         int width, int height,
                         const uint32_t *img, int img_w, int img_h,
                         float t)
{
    float scale = 0.016f;
    float half_w = img_w * scale * 0.5f;
    float half_h = img_h * scale * 0.5f;
    float depth_scale = 0.03f;  /* Subtle depth for 3D effect */

    /* Camera positioning - standing image */
    float model_size = (img_w > img_h ? img_w : img_h) * scale;
    float z_offset = -(model_size * 2.5f);

    b3d_init(pixels, depth, width, height, 60.0f);
    b3d_set_camera(0.0f, 0.0f, z_offset, 0.0f, 0.0f, 0.0f);
    b3d_clear();

    b3d_reset();
    b3d_rotate_y(t * 0.4f);  /* Rotate around Y axis like a spinning photo */

    for (int y = 0; y < img_h - 1; ++y) {
        for (int x = 0; x < img_w - 1; ++x) {
            uint32_t c00 = img[y * img_w + x];
            uint32_t c10 = img[y * img_w + (x + 1)];
            uint32_t c01 = img[(y + 1) * img_w + x];
            uint32_t c11 = img[(y + 1) * img_w + (x + 1)];

            float lum00 = ((c00 >> 16) & 0xff) * 0.2126f +
                          ((c00 >> 8) & 0xff) * 0.7152f +
                          (c00 & 0xff) * 0.0722f;
            float lum10 = ((c10 >> 16) & 0xff) * 0.2126f +
                          ((c10 >> 8) & 0xff) * 0.7152f +
                          (c10 & 0xff) * 0.0722f;
            float lum01 = ((c01 >> 16) & 0xff) * 0.2126f +
                          ((c01 >> 8) & 0xff) * 0.7152f +
                          (c01 & 0xff) * 0.0722f;
            float lum11 = ((c11 >> 16) & 0xff) * 0.2126f +
                          ((c11 >> 8) & 0xff) * 0.7152f +
                          (c11 & 0xff) * 0.0722f;

            /* Depth based on luminance (brighter = closer) */
            float d00 = (lum00 / 255.0f - 0.5f) * depth_scale;
            float d10 = (lum10 / 255.0f - 0.5f) * depth_scale;
            float d01 = (lum01 / 255.0f - 0.5f) * depth_scale;
            float d11 = (lum11 / 255.0f - 0.5f) * depth_scale;

            float fx = (x * scale) - half_w;
            float fy = half_h - (y * scale);  /* Flip Y so head is up */
            float fx1 = fx + scale;
            float fy1 = fy - scale;

            /* Standing image: X=horizontal, Y=vertical, Z=depth */
            b3d_triangle(fx,  fy,  d00,
                         fx1, fy,  d10,
                         fx1, fy1, d11, c00);
            b3d_triangle(fx,  fy,  d00,
                         fx1, fy1, d11,
                         fx,  fy1, d01, c11);
        }
    }
}

int main(int argc, char **argv)
{
    int img_w = 0, img_h = 0;
    uint32_t *img = decode_lena(&img_w, &img_h);
    if (!img) {
        fprintf(stderr, "Failed to decode Lena image\n");
        return 1;
    }

    int width = 1280;
    int height = 720;
    const char *snapshot = get_snapshot_path(argc, argv);

    uint32_t *pixels = malloc((size_t)width * (size_t)height * sizeof(pixels[0]));
    b3d_depth_t *depth = malloc((size_t)width * (size_t)height * sizeof(depth[0]));
    if (!pixels || !depth) {
        fprintf(stderr, "Allocation failed\n");
        free(img);
        return 1;
    }

    if (snapshot) {
        render_frame(pixels, depth, width, height, img, img_w, img_h, 1.2f);
        write_png(snapshot, pixels, width, height);
        free(pixels);
        free(depth);
        free(img);
        return 0;
    }

    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window *window = SDL_CreateWindow("Lena 3D", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, width, height, 0);
    SDL_Renderer *renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_PRESENTVSYNC);
    SDL_Texture *texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STREAMING, width, height);
    SDL_SetWindowTitle(window, "Lena heightfield (b3d)");

    int quit = 0;
    while (!quit) {
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT ||
               (event.type == SDL_KEYDOWN && event.key.keysym.scancode == SDL_SCANCODE_ESCAPE)) {
                quit = 1;
                break;
            }
        }
        if (quit) break;

        float t = SDL_GetTicks() * 0.001f;
        render_frame(pixels, depth, width, height, img, img_w, img_h, t);

        SDL_RenderClear(renderer);
        SDL_UpdateTexture(texture, NULL, pixels, width * sizeof(uint32_t));
        SDL_RenderCopy(renderer, texture, NULL, NULL);
        SDL_RenderPresent(renderer);
    }

    free(pixels);
    free(depth);
    free(img);
    SDL_DestroyTexture(texture);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}
