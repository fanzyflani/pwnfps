#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <errno.h>

#include <assert.h>
#include <math.h>

#include <SDL.h>

#include <mmintrin.h>
#include <xmmintrin.h>
#include <emmintrin.h>

#ifndef WIN32
#include <signal.h>
#endif

//define inline 
#define EPSILON 0.0000000000001f

#define REFLECT_BLUR 0.03f
#define PLAYER_BBOX 0.2f
#define REFLECT 2
#define POSTPROC_BLUR 1

#define DEF_SCALE 3
#define DEF_RWIDTH 320
#define DEF_RHEIGHT 200
#define DEF_WIDTH ((DEF_SCALE)*(DEF_RWIDTH))
#define DEF_HEIGHT ((DEF_SCALE)*(DEF_RHEIGHT))

#define COL_CEIL  _mm_setr_ps(30.0f, 30.0f, 0.0f, 0.0f)
#define COL_FLOOR _mm_setr_ps(1.0f, 1.0f, 1.0f, 0.0f)
#define COL_WALL  _mm_setr_ps(0.8f, 0.8f, 1.0f, 0.0f)

const float isqrt2 = 0.7071067811865476f;

float sec_current = 0.0f;

enum
{
	FXP = 0,
	FZP,
	FXN,
	FZN,

	FYP, FYN,
};

typedef union vec4_s
{
	float a[4];
	int32_t ai[4];
	struct { float x, y, z, w; } v;
	struct { float b, g, r, a; } c;
	struct { float s, t, r, q; } t;
	__m128 m;
	__m128i mi;
} vec4;

typedef union mat4_s
{
	vec4 a[4];
	struct { vec4 x, y, z, w; } v;
	struct { vec4 b, g, r, a; } c;
	struct { vec4 s, t, r, q; } t;
} mat4;

typedef enum part_typ_e
{
	P_INVAL = 0,

	P_SPHERE,

	P_CSG_UNION, // A + B
	P_CSG_DIFF, // A - B
	P_CSG_INTER, // A * B

	P_MAX
} part_typ;

typedef union part_s part;
union part_s
{
	part_typ typ;

	struct {
		part_typ typ;
		float r;
		float refl;
		vec4 pos;
		vec4 col;
	} sph;

	struct {
		part_typ typ;
		part *a, *b;
	} csg;
};

typedef struct portal_s
{
	int x1, z1;
	int x2, z2;
	int rot12;
	char c1, c2;
	int double1, double2;
} portal;

typedef struct level_s
{
	int sx, sz;

	portal pmap[26];

	char data[64][64];

	part **parts[64][64];
	uint16_t parts_num[64][64];
	uint16_t parts_max[64][64];

	part **objs;
	uint16_t objs_num;
} level;

int rwidth = DEF_RWIDTH;
int rheight = DEF_RHEIGHT;
int rscale = DEF_SCALE;

SDL_Surface *screen = NULL;
uint32_t *sbuf = NULL;
uint32_t *tsbuf = NULL;
float *zbuf = NULL;
level *lvroot = NULL;

static uint32_t randi(uint32_t *seed)
{
	*seed = (*seed * 25739) + 4;
	*seed &= 0x7FFFFFFF;
	return *seed;
}

static float randfu(uint32_t *seed)
{
	return (randi(seed)) % 3759 / 3759.0f;
}

static float randfs(uint32_t *seed)
{
	return (randfu(seed)) * 2.0f - 1.0f;
}

static inline float v_dot(__m128 va, __m128 vb)
{
	__m128 v2 = _mm_mul_ps(va, vb);

	// get that horizontal action working
	// TODO: get a better non-shufps method working
	__m128 va1 = _mm_add_ps(v2,
		_mm_shuffle_ps(v2, v2, 0x4E));
	__m128 va2 = _mm_add_ps(va1,
		_mm_shuffle_ps(va1, va1, 0xB1));

	return _mm_cvtss_f32(va2);
}

static inline __m128 v_normalise(__m128 v)
{
	__m128 v2 = _mm_mul_ps(v, v);

	// get that horizontal action working
	// TODO: get a better non-shufps method working
	__m128 va1 = _mm_add_ps(v2,
		_mm_shuffle_ps(v2, v2, 0x4E));
	__m128 va2 = _mm_add_ps(va1,
		_mm_shuffle_ps(va1, va1, 0xB1));

	__m128 vimag = _mm_rsqrt_ps(va2);

	return _mm_mul_ps(v, vimag);
}

uint32_t col_ftoint(__m128 v)
{
	__m128 vscale = _mm_mul_ps(v, _mm_set1_ps(255.0f));
	__m128i vint4 = _mm_cvtps_epi32(vscale);
	__m128i vint2 = _mm_packs_epi32(vint4, vint4);
	__m128i vint1 = _mm_packus_epi16(vint2, vint2);

	// TODO: find a better intrinsic for this
	uint32_t ret;
	_mm_store_ss((float *)&ret, (__m128)vint1);
	return ret;
}

void mat4_iden(mat4 *A)
{
	int i, j;

	for(i = 0; i < 4; i++)
	for(j = 0; j < 4; j++)
		A->a[i].a[j] = (i == j ? 1.0f : 0.0f);
}

void mat4_iden3(mat4 *A)
{
	int i, j;

	for(i = 0; i < 3; i++)
	for(j = 0; j < 3; j++)
		A->a[i].a[j] = (i == j ? 1.0f : 0.0f);
}


void mat4_rotx(mat4 *A, float ang)
{
	float vs = sinf(ang);
	float vc = cosf(ang);

	float vyy = A->v.y.v.y;
	float vyz = A->v.y.v.z;
	float vzy = A->v.z.v.y;
	float vzz = A->v.z.v.z;

	A->v.y.v.y = vc*vyy + vs*vyz;
	A->v.y.v.z = vc*vyz - vs*vyy;
	A->v.z.v.y = vc*vzy + vs*vzz;
	A->v.z.v.z = vc*vzz - vs*vzy;
}

void mat4_roty(mat4 *A, float ang)
{
	float vs = sinf(ang);
	float vc = cosf(ang);

	float vxx = A->v.x.v.x;
	float vxz = A->v.x.v.z;
	float vzx = A->v.z.v.x;
	float vzz = A->v.z.v.z;

	A->v.x.v.x = vc*vxx + vs*vxz;
	A->v.x.v.z = vc*vxz - vs*vxx;
	A->v.z.v.x = vc*vzx + vs*vzz;
	A->v.z.v.z = vc*vzz - vs*vzx;
}

static int celltype_is_solid(level *lv, char c, char oldcell, float y)
{
	if(c == '"' && (oldcell == '#' || oldcell == '&'))
		return y < 1.0f || y >= 2.0f;
	if(c == '#' || c == '&')
		return y < 0.0f || y >= 2.0f;
	if(c == ';' || c == '$' || c == '"')
		return y < 0.0f || y >= 1.0f;
	if(c == '>' || c == '<' || c == '^' || c == ',')
		return y < 0.0f || y >= 1.0f;
	if(c >= 'A' && c <= 'Z')
		return lv->pmap[c - 'A'].x2 != -1 ? 0 : 1;
	
	return 1;
}

static int celltype_is_free(char c)
{
	if(c == ';' || c == '$' || c == '"')
		return 1;
	if(c == '#' || c == '&')
		return 1;
	if(c == '>' || c == '<' || c == '^' || c == ',')
		return 1;
	
	return 0;
}

static int find_free_dir_2d(level *lv, int x, int z)
{
	if(celltype_is_free(lv->data[z][x+1])) return FXP;
	if(celltype_is_free(lv->data[z+1][x])) return FZP;
	if(celltype_is_free(lv->data[z][x-1])) return FXN;
	if(celltype_is_free(lv->data[z-1][x])) return FZN;

	printf("NOT FREE %i %i\n", x, z);
	return FXP; // stuff it
}

static char get_cell(level *lv, int cx, int cz)
{
	// Clamp for safety
	if(cx < 0 || cx >= 64) cx = 0;
	if(cz < 0 || cz >= 64) cz = 0;
	
	return lv->data[cz][cx];
}

static __m128 trace_ray(int hitctr, uint32_t *seed, level *lv, float *dist, const vec4 *ifrom, const vec4 *iray, __m128 icol);

static __m128 trace_hit_bounce(int hitctr, uint32_t *seed, level *lv, const vec4 *iray, const vec4 *ipos, const vec4 *inorm, int ldir, float refl, float *dist, float fog, __m128 col)
{
	//
	if(hitctr < 0 || hitctr >= REFLECT || refl == 0.0f)
		return col;

	vec4 ray, pos;
	ray.m = iray->m;
	pos.m = ipos->m;
	switch(ldir)
	{
		case FXP:
			ray.v.x = -ray.v.x;
			pos.v.x -= 0.001f;
			break;

		case FXN:
			ray.v.x = -ray.v.x;
			pos.v.x += 0.001f;
			break;

		case FZP:
			ray.v.z = -ray.v.z;
			pos.v.z -= 0.001f;
			break;

		case FZN:
			ray.v.z = -ray.v.z;
			pos.v.z += 0.001f;
			break;

		case FYP:
			ray.v.y = -ray.v.y;
			pos.v.y -= 0.001f;
			break;

		case FYN: {
			pos.v.y -= 0.001f;
			vec4 norm;
			float ang = ((float)M_PI)*2.0f*(
				sinf(((float)M_PI)*0.5f*pos.v.x)
				+ cosf(((float)M_PI)*0.5f*pos.v.z)
				+ sec_current*1.0f); 
			norm.m = v_normalise(_mm_setr_ps(sinf(ang), 38.0f, cosf(ang), 0.0f));

			float rmul = -2.0f * (0.0f
				+ ray.v.x*norm.v.x
				+ ray.v.y*norm.v.y
				+ ray.v.z*norm.v.z);

			ray.m = v_normalise(_mm_add_ps(
				_mm_mul_ps(_mm_set1_ps(rmul), norm.m),
				ray.m));
		} break;

		default: {
			pos.m = _mm_sub_ps(pos.m,
				_mm_mul_ps(_mm_set1_ps(0.001f), ray.m));

			vec4 norm;
			norm.m = inorm->m;

			float rmul = -2.0f * (0.0f
				+ ray.v.x*norm.v.x
				+ ray.v.y*norm.v.y
				+ ray.v.z*norm.v.z);

			ray.m = v_normalise(_mm_add_ps(
				_mm_mul_ps(_mm_set1_ps(rmul), norm.m),
				ray.m));
		} break;

	}

#ifdef REFLECT_BLUR
	const float rblur = REFLECT_BLUR;
	ray.v.x += randfs(seed) * rblur;
	ray.v.y += randfs(seed) * rblur;
	randfs(seed);
	ray.v.z += randfs(seed) * rblur;
	randfs(seed);
#endif

	float odist = *dist;
	*dist = 0;

	__m128 bcol = col;
	col = trace_ray(hitctr+1, seed, lv, dist, &pos, &ray, col);
	col = _mm_add_ps(
		_mm_mul_ps(_mm_set1_ps(refl), col),
		_mm_mul_ps(_mm_set1_ps(1.0f-refl), bcol));

	if(fog != 0.0f)
	{
		fog = expf(-0.6f*fog);
		col = _mm_add_ps(
			_mm_mul_ps(_mm_set1_ps(fog), col),
			_mm_mul_ps(_mm_set1_ps(1.0f-fog), _mm_set1_ps(1.0f)));
	}

	*dist = odist;

	return col;
}

static __m128 trace_hit_wall(int hitctr, uint32_t *seed, level *lv, const vec4 *ifrom, const vec4 *ipos, const vec4 *iray, int ldir, float *dist, __m128 icol, float fog, __m128 col)
{
	float diffuse;

	col = _mm_mul_ps(icol, col);

	switch(ldir)
	{
		default:
		case FXP:
			diffuse = iray->v.x;
			break;

		case FYP:
			diffuse = iray->v.y;
			break;

		case FZP:
			diffuse = iray->v.z;
			break;

		case FXN:
			diffuse = -iray->v.x;
			break;

		case FYN:
			diffuse = -iray->v.y;
			break;

		case FZN:
			diffuse = -iray->v.z;
			break;
	}

	if(diffuse < 0.0f)
		diffuse = 0.0f;

	// Apply ambient
	const float amb = 0.1f;
	diffuse = (1.0f - amb)*diffuse + amb;

	// Return colour
	col = _mm_mul_ps(col, _mm_set1_ps(diffuse));

	return trace_hit_bounce(hitctr, seed, lv, iray, ipos, NULL, ldir,
		(ldir == FYN ? 0.7f : 0.25f), dist, fog, col);
}

static void trace_ray_through(level *lv, int *ldir, float *cdist, vec4 *wdist, vec4 *pos, vec4 *ray, int gx, int gy, int gz)
{
	if(wdist->v.y < wdist->v.x && wdist->v.y < wdist->v.z)
	{
		// Y
		*cdist += wdist->v.y;
		pos->m = _mm_add_ps(_mm_mul_ps(
			_mm_set1_ps(wdist->v.y), ray->m),
			pos->m);
		*ldir = (gy < 0 ? FYN : FYP);

	} else if(wdist->v.x < wdist->v.z) {
		// X
		*cdist += wdist->v.x;
		pos->m = _mm_add_ps(_mm_mul_ps(
			_mm_set1_ps(wdist->v.x), ray->m),
			pos->m);
		*ldir = (gx < 0 ? FXN : FXP);

	} else {
		// Z
		*cdist += wdist->v.z;
		pos->m = _mm_add_ps(_mm_mul_ps(
			_mm_set1_ps(wdist->v.z), ray->m),
			pos->m);
		*ldir = (gz < 0 ? FZN : FZP);

	}
}

static __m128 trace_ray(int hitctr, uint32_t *seed, level *lv, float *dist, const vec4 *ifrom, const vec4 *iray, __m128 icol)
{
	int i;
	vec4 ray, pos;
	vec4 iavel, avel, wdist;
	int gx, gy, gz;
	float cdist = 0.0f;
	float fog = 0.0f;
	float fogbeg = 0.0f;

	float aux_dist = -1.0f;
	float aux_refl = 0.25f;
	int aux_dir = -1;
	vec4 aux_pos, aux_norm, aux_col;
	aux_col.m = _mm_set1_ps(1.0f);

	// Copy our constant inputs
	ray.m = iray->m;
	pos.m = ifrom->m;

	// Fuzz + normalise the direction
	/*
	ray.v.x += randfs(seed) * 0.004f;
	ray.v.y += randfs(seed) * 0.004f;
	ray.v.z += randfs(seed) * 0.004f;
	*/
	ray.m = v_normalise(ray.m);

	// Find our starting point
	int cx = (int)ifrom->v.x;
	int cz = (int)ifrom->v.z;

	// Move velocities away from zero
	// TODO: SSE version
	if(ray.v.x > -EPSILON && ray.v.x < EPSILON) ray.v.x = (ray.v.x < 0.0f ? -EPSILON : EPSILON);
	if(ray.v.y > -EPSILON && ray.v.y < EPSILON) ray.v.y = (ray.v.y < 0.0f ? -EPSILON : EPSILON);
	if(ray.v.z > -EPSILON && ray.v.z < EPSILON) ray.v.z = (ray.v.z < 0.0f ? -EPSILON : EPSILON);

	// Get direction
	gx = (iray->v.x < 0 ? -1 : 1);
	gy = (iray->v.y < 0 ? -1 : 1);
	gz = (iray->v.z < 0 ? -1 : 1);

	// Get abs velocity
	avel.m = (__m128)_mm_and_si128(_mm_set1_epi32(0x7FFFFFFF), (__m128i)ray.m);
	iavel.m = _mm_rcp_ps(avel.m);

	// Get distance to cell boundary
	// TODO: proper SSE version
	wdist.m = _mm_sub_ps(pos.m, _mm_cvtepi32_ps(_mm_setr_epi32(cx, 0, cz, 0)));
	if(ray.v.x >= 0.0f) wdist.v.x = 1.0f - wdist.v.x; 
	if(ray.v.y >= 0.0f) wdist.v.y = 1.0f - wdist.v.y; 
	if(ray.v.z >= 0.0f) wdist.v.z = 1.0f - wdist.v.z; 

	// Now convert distance to time
	wdist.m = _mm_mul_ps(wdist.m, iavel.m);

	// Get a pointer + the level width
	char cell = get_cell(lv, cx, cz);

	// Do our trace!
	int maxsteps = 1000;
	int ldir = FYN;

	while(maxsteps-- > 0)
	{
		if(cx >= 0 && cx < 64 && cz >= 0 && cz < 64)
		for(i = 0; i < lv->parts_num[cz][cx]; i++)
		{
			part *pt = lv->parts[cz][cx][i];

			// TODO: handle stuff other than spheres

			/*
			pt->typ = P_SPHERE;
			pt->sph.r = 0.3f;
			//sph_rpos.m = _mm_setr_ps(8.5f, 0.3f, 7.0f, 1.0f);
			//sph_rpos.m = _mm_setr_ps(fmodf(0.5f + sec_current/2.0f, 1.0f), 0.3f, 0.5f, 1.0f);
			pt->sph.pos.m = _mm_setr_ps(0.5f, 0.8f, 0.5f, 1.0f);
			pt->sph.pos.v.x += (float)cx;
			pt->sph.pos.v.z += (float)cz;
			pt->sph.col.m = _mm_setr_ps(0.8f, 1.0f, 0.8f, 1.0f);
			*/

			float sph_rad2 = pt->sph.r * pt->sph.r;
			float sph_dist = -1.0f;

			// Get relative pos
			vec4 sph_pos;
			sph_pos.m = _mm_sub_ps(pt->sph.pos.m, pos.m);

			// Get distance + dot
			float sph_dist2 = v_dot(sph_pos.m, sph_pos.m);
			float sph_dot = v_dot(sph_pos.m, ray.m);
			if(sph_dot > 0.0f)
			{
				float sph_calcrad2 = (sph_dist2 - sph_dot*sph_dot);

				if(sph_calcrad2 < sph_rad2)
				{
					// Get correct dist
					float sph_sdist2 = 1.0f - sph_calcrad2/sph_rad2;
					sph_dist = sqrtf(sph_dist2) - sqrtf(sph_sdist2);

					if(aux_dist == -1.0f || sph_dist + cdist < aux_dist)
					{
						aux_dist = sph_dist + cdist;
						aux_pos.m = _mm_add_ps(pos.m,
							_mm_mul_ps(_mm_set1_ps(sph_dist), ray.m));

						// Do diffuse
						aux_norm.m = v_normalise(_mm_sub_ps(aux_pos.m, pt->sph.pos.m));
						float diff = -v_dot(ray.m, aux_norm.m);
						if(diff < 0.0f) diff = 0.0f;
						float amb = 0.2f;
						aux_refl = pt->sph.refl;
						diff = amb + (1.0f - amb)*diff;
						aux_col.m = _mm_mul_ps(_mm_set1_ps(diff), pt->sph.col.m);
					}
				}
			}
		}

		// Work out what to do with this cell
		char this_cell = cell;
		switch(this_cell)
		{
			case ';':
			case '$':
			case '"':
				// Regular box.
				// Trace to edge

				//printf("%f %f %f %f\n", wdist.v.x, wdist.v.z, iavel.v.x, iavel.v.z);
				if(this_cell == '$') fogbeg = cdist;

				trace_ray_through(lv, &ldir, &cdist, &wdist, &pos, &ray, gx, gy, gz);

				if(aux_dist != -1.0f && cdist > aux_dist)
				{
					*dist = aux_dist;
					if(this_cell == '$' && aux_dist > fogbeg) fog += aux_dist - fogbeg;
					return trace_hit_bounce(hitctr, seed, lv, &ray, &aux_pos, &aux_norm, aux_dir, aux_refl,
						dist, fog, aux_col.m);
				}

				if(this_cell == '$') fog += cdist - fogbeg;

				if(ldir == FYN || ldir == FYP)
				{
					*dist = cdist;
					return trace_hit_wall(hitctr, seed, lv, ifrom, &pos, &ray, ldir, dist, icol, fog,
						gy > 0
						? COL_CEIL
						: COL_FLOOR);

				} else if(ldir == FXN || ldir == FXP) {
					wdist.m = _mm_sub_ps(wdist.m, _mm_set1_ps(wdist.v.x));
					wdist.v.x = iavel.v.x;
					cx += gx;

				} else if(ldir == FZN || ldir == FZP) {
					wdist.m = _mm_sub_ps(wdist.m, _mm_set1_ps(wdist.v.z));
					wdist.v.z = iavel.v.z;
					cz += gz;
				}

				cell = get_cell(lv, cx, cz);
				if(this_cell == '"' && (cell == '#' || cell == '&'))
				{
					pos.v.y += 1.0f;

					if(gy < 0.0f)
						wdist.v.y += iavel.v.y;
					else
						wdist.v.y -= iavel.v.y;
				}
				break;

			case '#':
			case '&':
				// 2-high room
				if(gy > 0) wdist.v.y += iavel.v.y;

				if(this_cell == '&') fogbeg = cdist;

				trace_ray_through(lv, &ldir, &cdist, &wdist, &pos, &ray, gx, gy, gz);

				if(aux_dist != -1.0f && cdist > aux_dist)
				{
					*dist = aux_dist;
					if(this_cell == '&' && aux_dist > fogbeg) fog += aux_dist - fogbeg;
					return trace_hit_bounce(hitctr, seed, lv, &ray, &aux_pos, &aux_norm, aux_dir, aux_refl,
						dist, fog, aux_col.m);
				}

				if(this_cell == '&') fog += cdist - fogbeg;

				if(ldir == FYN || ldir == FYP)
				{
					*dist = cdist;
					return trace_hit_wall(hitctr, seed, lv, ifrom, &pos, &ray, ldir, dist, icol, fog,
						gy > 0
						? COL_CEIL
						: COL_FLOOR);

				} else if(ldir == FXN || ldir == FXP) {
					wdist.m = _mm_sub_ps(wdist.m, _mm_set1_ps(wdist.v.x));
					wdist.v.x = iavel.v.x;
					cx += gx;

				} else if(ldir == FZN || ldir == FZP) {
					wdist.m = _mm_sub_ps(wdist.m, _mm_set1_ps(wdist.v.z));
					wdist.v.z = iavel.v.z;
					cz += gz;
				}

				if(gy > 0) wdist.v.y -= iavel.v.y;

				cell = get_cell(lv, cx, cz);
				if(cell == '"')
				{
					pos.v.y -= 1.0f;
					if(gy > 0.0f)
						wdist.v.y += iavel.v.y;
					else
						wdist.v.y -= iavel.v.y;
				}

				char xcell = cell;
				if(xcell >= 'A' && xcell <= 'Z')
				{
					portal *pm = &lv->pmap[xcell - 'A'];

					if(pm->x1 == cx && pm->z1 == cz)
						xcell = pm->c2;
					else if(pm->x2 == cx && pm->z2 == cz)
						xcell = pm->c1;
				}

				if(pos.v.y < 0.0f || pos.v.y > 1.0f)
				switch(xcell)
				{
					case '#':
					case '&':
						break;

					case '"':
						// revert our adjustments
						// otherwise reflections are broken
						pos.v.y += 1.0f;

						if(gy > 0.0f)
							wdist.v.y -= iavel.v.y;
						else
							wdist.v.y += iavel.v.y;

						*dist = cdist;
						return trace_hit_wall(hitctr, seed, lv, ifrom, &pos, &ray, ldir, dist, icol, fog,
							COL_WALL);

					default:
						*dist = cdist;
						return trace_hit_wall(hitctr, seed, lv, ifrom, &pos, &ray, ldir, dist, icol, fog,
							COL_WALL);
				}
				break;

			case '>':
			case '<':
			case ',':
			case '^': {
				// Ramp up
				// TODO: ? make this not distort from the sides?
				// Tilt ray
				const float ramp_delta = 0.5f;
				switch(this_cell)
				{
					case '>': ray.v.y -= ramp_delta * ray.v.x; break;
					case '<': ray.v.y += ramp_delta * ray.v.x; break;
					case ',': ray.v.y -= ramp_delta * ray.v.z; break;
					case '^': ray.v.y += ramp_delta * ray.v.z; break;
				}

				wdist.v.y = pos.v.y;
				if(ray.v.y >= 0.0f) wdist.v.y = 1.0f - wdist.v.y; 
				wdist.v.y *= 1.0f / (ray.v.y < 0.0f ? -ray.v.y : ray.v.y);

				if(aux_dist != -1.0f && cdist > aux_dist)
				{
					*dist = aux_dist;
					return trace_hit_bounce(hitctr, seed, lv, &ray, &aux_pos, &aux_norm, aux_dir, aux_refl,
						dist, fog, aux_col.m);
				}

				trace_ray_through(lv, &ldir, &cdist, &wdist, &pos, &ray, gy, gy, gz);
				if(ldir == FYN || ldir == FYP)
				{
					*dist = cdist;
					ldir = (ray.v.y < 0.0f ? FYN : FYP);
					return trace_hit_wall(hitctr, seed, lv, ifrom, &pos, &ray, ldir, dist, icol, fog,
						ray.v.y >= 0.0f
						? COL_CEIL
						: COL_FLOOR);

				} else if(ldir == FXN || ldir == FXP) {
					ldir = (ray.v.x < 0.0f ? FXN : FXP);
					wdist.m = _mm_sub_ps(wdist.m, _mm_set1_ps(wdist.v.x));
					wdist.v.x = iavel.v.x;
					cx += gx;

				} else if(ldir == FZN || ldir == FZP) {
					ldir = (ray.v.z < 0.0f ? FZN : FZP);
					wdist.m = _mm_sub_ps(wdist.m, _mm_set1_ps(wdist.v.z));
					wdist.v.z = iavel.v.z;
					cz += gz;
				}

				switch(this_cell)
				{
					case '>': ray.v.y += ramp_delta * ray.v.x; break;
					case '<': ray.v.y -= ramp_delta * ray.v.x; break;
					case ',': ray.v.y += ramp_delta * ray.v.z; break;
					case '^': ray.v.y -= ramp_delta * ray.v.z; break;
				}

				wdist.v.y = pos.v.y;
				if(ray.v.y >= 0.0f) wdist.v.y = 1.0f - wdist.v.y; 
				wdist.v.y *= iavel.v.y;
				cell = get_cell(lv, cx, cz);
			} break;

			default: // Treat this as a wall.
			if(cell >= 'A' && cell <= 'Z')
			{
				portal *pm = &lv->pmap[cell - 'A'];

				if(pm->x2 == -1)
				{
					// ERROR
					*dist = cdist;
					if(aux_dist != -1.0f && cdist > aux_dist)
					{
						*dist = aux_dist;
						return trace_hit_bounce(hitctr, seed, lv, &ray, &aux_pos, &aux_norm, aux_dir, aux_refl,
							dist, fog, aux_col.m);
					}

					return trace_hit_wall(hitctr, seed, lv, ifrom, &pos, &ray, ldir, dist, icol, fog,
						//_mm_setr_ps(0.0f, 0.0f, 5.0f, 0.0f));
						COL_WALL);
				}

				// Find out which cell we are
				int rot = 0;
				if(pm->x1 == cx && pm->z1 == cz)
				{
					// 1 -> 2
					cx = pm->x2;
					cz = pm->z2;
					pos.v.x += pm->x2 - pm->x1;
					pos.v.z += pm->z2 - pm->z1;
					rot = (-pm->rot12) & 3;

				} else if(pm->x2 == cx && pm->z2 == cz) {
					// 2 -> 1
					cx = pm->x1;
					cz = pm->z1;
					pos.v.x -= pm->x2 - pm->x1;
					pos.v.z -= pm->z2 - pm->z1;
					rot = (pm->rot12) & 3;

				} else {
					// ERROR
					*dist = cdist;
					if(aux_dist != -1.0f && cdist > aux_dist)
					{
						*dist = aux_dist;
						return trace_hit_bounce(hitctr, seed, lv, &ray, &aux_pos, &aux_norm, aux_dir, aux_refl,
							dist, fog, aux_col.m);
					}

					return trace_hit_wall(hitctr, seed, lv, ifrom, &pos, &ray, ldir, dist, icol, fog,
						_mm_setr_ps(5.0f, 0.0f, 5.0f, 0.0f));
				}

				float rx = pos.v.x;
				float rz = pos.v.z;
				int tgx = gx;
				int tgz = gz;
				float rvzx = ray.v.x;
				float rvzz = ray.v.z;
				float trx = rx;
				float trz = rz;
				float trvzx = rvzx;
				float trvzz = rvzz;
				float twdx = wdist.v.x;
				float twdz = wdist.v.z;
				float tiax = iavel.v.x;
				float tiaz = iavel.v.z;

				ldir = (ldir - rot) & 3;

				switch(rot)
				{
					case 0:
						break;

					case 1:
						rx = (cx + 0.5f) + (trz-(cz + 0.5f));
						rz = (cz + 0.5f) - (trx-(cx + 0.5f));
						rvzx =  trvzz;
						rvzz = -trvzx;
						gx =  tgz;
						gz = -tgx;
						wdist.v.x = twdz;
						wdist.v.z = twdx;
						iavel.v.x = tiaz;
						iavel.v.z = tiax;
						break;

					case 2:
						rx = (cx + 0.5f)*2.0f - rx;
						rz = (cz + 0.5f)*2.0f - rz;
						rvzx = -trvzx;
						rvzz = -trvzz;
						gx = -gx;
						gz = -gz;
						break;

					case 3:
						rx = (cx + 0.5f) - (trz-(cz + 0.5f));
						rz = (cz + 0.5f) + (trx-(cx + 0.5f));
						rvzx = -trvzz;
						rvzz =  trvzx;
						gx = -tgz;
						gz =  tgx;
						wdist.v.x = twdz;
						wdist.v.z = twdx;
						iavel.v.x = tiaz;
						iavel.v.z = tiax;
						break;
				}

				pos.v.x = rx;
				pos.v.z = rz;
				ray.v.x = rvzx;
				ray.v.z = rvzz;

				switch(ldir)
				{
					default:

					case FXP:
						cx++;
						pos.v.x += 1.0f;
						break;

					case FZP:
						cz++;
						pos.v.z += 1.0f;
						break;

					case FXN:
						cx--;
						pos.v.x -= 1.0f;
						break;

					case FZN:
						cz--;
						pos.v.z -= 1.0f;
						break;
				}

				cell = get_cell(lv, cx, cz);
				break;
			} else {
				*dist = cdist;
				if(aux_dist != -1.0f && cdist > aux_dist)
				{
					*dist = aux_dist;
					return trace_hit_bounce(hitctr, seed, lv, &ray, &aux_pos, &aux_norm, aux_dir, aux_refl,
						dist, fog, aux_col.m);
				}

				return trace_hit_wall(hitctr, seed, lv, ifrom, &pos, &ray, ldir, dist, icol, fog,
					ldir == FYP
					? COL_CEIL
					: COL_WALL);
			}

		}

		if(aux_dist != -1.0f && cdist > aux_dist)
		{
			*dist = aux_dist;
			return trace_hit_bounce(hitctr, seed, lv, &ray, &aux_pos, &aux_norm, aux_dir, aux_refl,
				dist, fog, aux_col.m);
		}

	}

	// TODO: handle OOB properly
	return ray.m;
}

// [p1,p2)
void trace_screen_centred(level *lv, int x1, int y1, int x2, int y2, const mat4 *cam)
{
	assert(x1 < x2);
	assert(y1 < y2);

	// Calculate dimensions
	int dimx = x2 - x1;
	int dimy = y2 - y1;

	// Calculate ratios
	float xrat = -1.0f;
	float yrat = -((float)dimy)/(float)dimx;
	float xsrat = 2.0f * xrat / (float)dimx;
	float ysrat = 2.0f * yrat / (float)dimy;

	// Work out x, y stuff
	vec4 rayb, rdx, rdy;
	const vec4 *from = &cam->v.w;

	rayb.m = _mm_add_ps(cam->v.z.m,
		_mm_add_ps(
			_mm_mul_ps(_mm_set1_ps(-xrat), cam->v.x.m),
			_mm_mul_ps(_mm_set1_ps(-yrat), cam->v.y.m)));

	rdx.m = _mm_mul_ps(_mm_set1_ps(xsrat), cam->v.x.m);
	rdy.m = _mm_mul_ps(_mm_set1_ps(ysrat), cam->v.y.m);

	// Trace!
	// Do 32x32 tiles for improved cache coherence
	int cy;
#pragma omp parallel for
	for(cy = 0; cy < dimy; cy += 32)
	{
		int y, cx, x;

		for(cx = 0; cx < dimx; cx += 32)
		for(y = 0; y < 32 && y + cy < dimy; y++)
		{
			vec4 rayl;
			rayl.m = _mm_add_ps(rayb.m, _mm_mul_ps(_mm_set1_ps(y+cy), rdy.m));
			rayl.m = _mm_add_ps(rayl.m, _mm_mul_ps(_mm_set1_ps(cx), rdx.m));
			uint32_t *p = sbuf + (y+cy)*rwidth + cx;
			float *dist = zbuf + (y+cy)*rwidth + cx;
			for(x = 0; x < 32 && x + cx < dimx; x++)
			{
				rayl.m = _mm_add_ps(rayl.m, rdx.m);
				uint32_t seed = (x+cx) + (y+cy)*(y+cy)*(rwidth+1);
				seed *= seed * seed;
				seed *= seed * seed;
				*(p++) = col_ftoint(
					trace_ray(0, &seed, lv, dist++, from, &rayl,
					_mm_setr_ps(1.0f, 1.0f, 1.0f, 1.0f)));
			}
		}
	}

#ifdef POSTPROC_BLUR
	int blurctr;
	for(blurctr = 0; blurctr < POSTPROC_BLUR; blurctr++)
	{
		// Apply focal blur
		// TODO: isolate this to the camera region
		memcpy(tsbuf, sbuf, sizeof(uint32_t)*rwidth*rheight);

	#pragma omp parallel for
		for(cy = 0; cy < dimy; cy++)
		{
			int cx, x, y;
			int i;
			uint32_t seed = cy*cy + 415135; // keyboard mash
			float *dist = zbuf + cy*rwidth;
			uint32_t *p = sbuf + cy*rwidth;

			const float fstr = 0.002f * dimy;
			const float foffs = 1.0f;
			// XXX: WARNING: EXPECTS A WIDTH DIVISIBLE BY 4
			// i'm currently just not blurring the last dimx%4 pixels in that case

			for(cx = 0; cx < dimx-3; cx += 4)
			{
				vec4 vacc;
				vec4 vbuf[4];
				for(i = 0; i < 4; i++)
				{
					int j;
					for(j = 0; j < 4; j++)
					{
						float z = (dist[j] - foffs);
						x = cx + j + randfs(&seed) * fstr * z;
						y = cy + randfs(&seed) * fstr * z;
						if(x < 0) x = 0;
						if(y < 0) y = 0;
						if(x >= dimx) x = dimx-1;
						if(y >= dimy) y = dimy-1;
						vbuf[i].ai[j] = tsbuf[x + y*rwidth];
					}
				}


				vacc.mi = _mm_avg_epu8(
					_mm_avg_epu8(vbuf[0].mi, vbuf[1].mi),
					_mm_avg_epu8(vbuf[2].mi, vbuf[3].mi)
				);

				_mm_store_si128((void *)p, vacc.mi);
				dist += 4;
				p += 4;
			}
		}
	}
#endif
}

void level_part_add_bbox(level *lv, part *pt, int cx1, int cz1, int cx2, int cz2)
{
	int x, z;

	// TODO: handle portal transitions
	for(z = cz1; z <= cz2; z++)
	for(x = cx1; x <= cx2; x++)
	{
		int idx = lv->parts_num[z][x]++;
		if(idx >= lv->parts_max[z][x])
		{
			lv->parts_max[z][x] = idx + 3;
			lv->parts[z][x] = realloc(lv->parts[z][x], lv->parts_max[z][x]*sizeof(part *));
		}

		lv->parts[z][x][idx] = pt;
	}

}

void level_part_add(level *lv, part *pt)
{
	switch(pt->typ)
	{
		case P_SPHERE:
			level_part_add_bbox(lv, pt,
				pt->sph.pos.v.x - pt->sph.r,
				pt->sph.pos.v.z - pt->sph.r,
				pt->sph.pos.v.x + pt->sph.r,
				pt->sph.pos.v.z + pt->sph.r);
			break;

		default:
			printf("unsupported part type %i\n", pt->typ);
			fflush(stdout);
			abort();
	}

}

void level_prepare_render(level *lv)
{
	int x, z;

	// Clear parts grid
	for(z = 0; z < 64; z++)
	for(x = 0; x < 64; x++)
		lv->parts_num[z][x] = 0;
}

void screen_upscale(void)
{
	int px, py, x, y;
	int pitch = screen->pitch/sizeof(uint32_t);
	int w = rwidth;
	int h = rheight;
	int ypitch = pitch*(rscale-1);
	int fypitch = pitch*rscale;

	uint32_t *s = sbuf;
	uint32_t *d = (uint32_t *)screen->pixels;

	for(py = 0; py < h; py++, d += ypitch)
	for(px = 0; px < w; px++, d += rscale)
	{
		uint32_t pd = *(s++);

		for(y = 0; y < rscale; y++, d += pitch)
		for(x = 0; x < rscale; x++)
			d[x] = pd;

		d -= fypitch;
	}
}

level *level_load(const char *fname)
{
	int x, z, i;
	FILE *fp = fopen(fname, "rb");
	if(fp == NULL)
	{
		printf("err loading level: %s\n", strerror(errno));
		return NULL;
	}

	char *tdata = malloc(64*64);
	level *lv = malloc(sizeof(level));
	memset(tdata, '.', 64*64);
	lv->objs = NULL;
	lv->objs_num = 0;

	for(i = 0; i < 26; i++)
	{
		lv->pmap[i].x1 = lv->pmap[i].x2 = -1;
		lv->pmap[i].c1 = ';';
		lv->pmap[i].c2 = ';';
	}

	lv->sx = 0;
	lv->sz = 0;
	for(z = 0; z < 64; z++)
	{
		for(x = 0; x < 64; x++)
		{
			int c = fgetc(fp);

			if(c == -1 || c == '\r' || c == '\n')
			{
				if(c == -1)
				{
					goto done_load;
				} else if(x == 0) {
					x--;
					continue;
				} else {
					break;
				}
			}

			if(c == '*')
			{
				c = ';';
				lv->sx = x;
				lv->sz = z;
			}

			if(c >= 'a' && c <= 'z'-1)
			{
				//
				portal *pm = &lv->pmap[c-'a'];

				if(pm->x1 == -1)
				{
					pm->x1 = x;
					pm->z1 = z;

				} else if(pm->x2 == -1) {
					pm->x2 = x;
					pm->z2 = z;

				}

				c = (c-'a')+'A' + 1;
			}

			if(c >= 'A' && c <= 'Z')
			{
				//
				portal *pm = &lv->pmap[c-'A'];

				if(pm->x1 == -1)
				{
					pm->x1 = x;
					pm->z1 = z;

				} else if(pm->x2 == -1) {
					pm->x2 = x;
					pm->z2 = z;

				}
			}

			tdata[x + z*64] = c;
		}
	}

	done_load:

	for(z = 0; z < 64; z++)
	for(x = 0; x < 64; x++)
	{
		lv->data[z][x] = tdata[x + z*64];
		lv->parts[z][x] = NULL;
		lv->parts_num[z][x] = 0;
		lv->parts_max[z][x] = 0;
	}
	
	for(i = 0; i < 26; i++)
	{
		portal *pm = &lv->pmap[i];

		if(pm->x2 == -1) continue;

		int d1 = find_free_dir_2d(lv, pm->x1, pm->z1);
		int d2 = find_free_dir_2d(lv, pm->x2, pm->z2);
		pm->rot12 = (d2 - d1 + 2) & 3;
		
		switch(d1)
		{
			case FXP: pm->c1 = lv->data[pm->z1][pm->x1+1]; break;
			case FZP: pm->c1 = lv->data[pm->z1+1][pm->x1]; break;
			case FXN: pm->c1 = lv->data[pm->z1][pm->x1-1]; break;
			case FZN: pm->c1 = lv->data[pm->z1-1][pm->x1]; break;
		}
		
		switch(d2)
		{
			case FXP: pm->c2 = lv->data[pm->z2][pm->x2+1]; break;
			case FZP: pm->c2 = lv->data[pm->z2+1][pm->x2]; break;
			case FXN: pm->c2 = lv->data[pm->z2][pm->x2-1]; break;
			case FZN: pm->c2 = lv->data[pm->z2-1][pm->x2]; break;
		}

		printf("%c rot %i\n", i+'A', pm->rot12);
	}

	free(tdata);

	fclose(fp);

	printf("spawn %3i %3i\n", lv->sx, lv->sz);

	return lv;
}

int mainloop(void)
{
	int i;

	SDL_Event ev;
	mat4 cam;
	vec4 gravity;

	lvroot = level_load("level.txt");
	level *lv = lvroot;

	gravity.m = _mm_setzero_ps();
	mat4_iden(&cam);
	cam.v.w.v.x = 0.5f + lv->sx;
	cam.v.w.v.y = 0.5f;
	cam.v.w.v.z = 0.5f + lv->sz;

	int tlast = (int)SDL_GetTicks();
	int tlastfps = tlast;
	int fps = 0;

	int k_turnleft = 0;
	int k_turnright = 0;
	int k_turnup = 0;
	int k_turndown = 0;
	int k_moveforward = 0;
	int k_moveback = 0;
	int k_moveleft = 0;
	int k_moveright = 0;

	float cangx = 0.0f;
	float cangy = 0.0f;

	part sph[3];
	for(i = 0; i < 3; i++)
	{
		sph[i].typ = P_SPHERE;
		sph[i].sph.r = 0.3f;
		sph[i].sph.refl = 0.6f;
		sph[i].sph.col.m = _mm_setr_ps(1.0f, 0.5f, 0.5f, 1.0f);
	}

	for(;;)
	{
		for(i = 0; i < 3; i++)
		{
			sph[i].sph.pos.m = _mm_setr_ps(
				2.0f+lv->sx + sinf(M_PI*(sec_current*0.2f + i*2/3.0f))*0.7f,
				0.3f + 0.1f + sinf(2.0f*M_PI*(sec_current*0.2f + i*2/3.0f))*0.1f,
				2.0f+lv->sz + cosf(M_PI*(sec_current*0.2f + i*2/3.0f))*0.7f,
				1.0f);
		}

		level_prepare_render(lv);
		level_part_add(lv, &sph[0]);
		level_part_add(lv, &sph[1]);
		level_part_add(lv, &sph[2]);
		trace_screen_centred(lv, 0, 0, rwidth, rheight, &cam);
		screen_upscale();
		SDL_Flip(screen);
		int tnext = (int)SDL_GetTicks();
		//SDL_Delay(10);

		float tdiff = (tnext - tlast) * (1.0f/1000.0f);
		sec_current += tdiff;
		tlast = tnext;
		fps++;

		if(tnext - tlastfps >= 1000)
		{
			tlastfps += 1000;
			char nbuf[128];
			sprintf(nbuf, "7DFPS 2014 - FPS: %i", fps);
			SDL_WM_SetCaption(nbuf, NULL);
			fps = 0;
		}

		while(SDL_PollEvent(&ev))
		switch(ev.type)
		{
			case SDL_QUIT:
				return 0;
			case SDL_KEYDOWN:
			case SDL_KEYUP:
			switch(ev.key.keysym.sym)
			{
				default:
					break;

				case SDLK_LEFT:
					k_turnleft = (ev.type == SDL_KEYDOWN ? 1 : 0);
					break;

				case SDLK_RIGHT:
					k_turnright = (ev.type == SDL_KEYDOWN ? 1 : 0);
					break;
				
				case SDLK_UP:
					k_turnup = (ev.type == SDL_KEYDOWN ? 1 : 0);
					break;

				case SDLK_DOWN:
					k_turndown = (ev.type == SDL_KEYDOWN ? 1 : 0);
					break;

				case SDLK_w:
					k_moveforward = (ev.type == SDL_KEYDOWN ? 1 : 0);
					break;

				case SDLK_s:
					k_moveback = (ev.type == SDL_KEYDOWN ? 1 : 0);
					break;

				case SDLK_a:
					k_moveleft = (ev.type == SDL_KEYDOWN ? 1 : 0);
					break;

				case SDLK_d:
					k_moveright = (ev.type == SDL_KEYDOWN ? 1 : 0);
					break;
			} break;
		}

		mat4_roty(&cam, tdiff*3.0f*(k_turnleft - k_turnright));
		//cangy += tdiff*3.0f*(k_turnleft - k_turnright);
		//cangx += tdiff*3.0f*(k_turnup - k_turndown);
		//mat4_iden3(&cam);
		//mat4_rotx(&cam, cangx); // TODO!
		//mat4_roty(&cam, cangy);

		// Take old cell + pos
		int cx1 = (int)cam.v.w.v.x;
		int cz1 = (int)cam.v.w.v.z;

		// Get velocity
		vec4 vel;
		vel.m = _mm_add_ps(
			_mm_mul_ps(cam.v.z.m,
				_mm_set1_ps(tdiff*5.0f*(k_moveforward - k_moveback))),
			_mm_mul_ps(cam.v.x.m,
				_mm_set1_ps(tdiff*5.0f*(k_moveleft - k_moveright))));

		// Move
		cam.v.w.m = _mm_add_ps(cam.v.w.m, vel.m);

		// Get newer pos
		float px1 = cam.v.w.v.x;
		float py1 = cam.v.w.v.y;
		float pz1 = cam.v.w.v.z;

		// Perform pushbacks
		int gx1 = (vel.v.x < 0.0f ? -1 : 1);
		int gy1 = (vel.v.y < 0.0f ? -1 : 1);
		int gz1 = (vel.v.z < 0.0f ? -1 : 1);
		float bx1 = px1 + gx1*PLAYER_BBOX;
		float by1 = py1 + gy1*PLAYER_BBOX;
		float bz1 = pz1 + gz1*PLAYER_BBOX;

		int bcx = (int)bx1;
		int bcy = (int)by1;
		int bcz = (int)bz1;
		char oldcell = get_cell(lv, cx1, cz1);

		if(cx1 != bcx && cz1 != bcz)
		{
			char cellx = get_cell(lv, bcx, cz1);
			char cellz = get_cell(lv, cx1, bcz);
			char cellc = get_cell(lv, bcx, bcz);

			int solx = celltype_is_solid(lv, cellx, oldcell, py1);
			int solz = celltype_is_solid(lv, cellz, oldcell, py1);
			int solc = celltype_is_solid(lv, cellc, oldcell, py1);

			if(solx && solz)
			{
				cam.v.w.v.x = cx1 + 0.5f + (0.5f-PLAYER_BBOX)*gx1;
				cam.v.w.v.z = cz1 + 0.5f + (0.5f-PLAYER_BBOX)*gz1;

			} else if(solx) {
				cam.v.w.v.x = cx1 + 0.5f + (0.5f-PLAYER_BBOX)*gx1;

			} else if(solz) {
				cam.v.w.v.z = cz1 + 0.5f + (0.5f-PLAYER_BBOX)*gz1;

			} else if(solc) {
				// TODO: quick trace to the wall
				cam.v.w.v.z = cz1 + 0.5f + (0.5f-PLAYER_BBOX)*gz1;

			}
		} else if(cx1 != bcx) {
			char c = get_cell(lv, bcx, bcz);

			if(celltype_is_solid(lv, c, oldcell, py1))
				cam.v.w.v.x = cx1 + 0.5f + (0.5f-PLAYER_BBOX)*gx1;
		} else if(cz1 != bcz) {
			char c = get_cell(lv, bcx, bcz);

			if(celltype_is_solid(lv, c, oldcell, py1))
				cam.v.w.v.z = cz1 + 0.5f + (0.5f-PLAYER_BBOX)*gz1;
		}

		// Apply gravity
		cam.v.w.m = _mm_add_ps(cam.v.w.m, gravity.m);
		gravity.v.y -= 3.0f*tdiff*tdiff;

		if(cam.v.w.v.y < 0.4f)
		{
			cam.v.w.v.y = 0.4f;
			gravity.v.y = 0.0f;
		}

		// Get new cell
		int cx2 = (int)cam.v.w.v.x;
		int cz2 = (int)cam.v.w.v.z;

		if(cx1 != cx2 || cz1 != cz2)
		{
			// Check cell types
			char c1 = get_cell(lv, cx1, cz1);
			char c2 = get_cell(lv, cx2, cz2);

			if((c1 == '#' || c1 == '&') && c2 == '"')
			{
				cam.v.w.v.y -= 1.0f;

			} else if(c1 == '"' && (c2 == '#' || c2 == '&')) {
				cam.v.w.v.y += 1.0f;

			} else if(c2 >= 'A' && c2 <= 'Z') {
				portal *pm = &lv->pmap[c2 - 'A'];

				printf("TRAVERSE %c\n", c2);

				int rot = 0;
				float rx = cam.v.w.v.x;
				float rz = cam.v.w.v.z;
				float rcx = cx2;
				float rcz = cz2;
				float rvxx = cam.v.x.v.x;
				float rvxz = cam.v.x.v.z;
				float rvzx = cam.v.z.v.x;
				float rvzz = cam.v.z.v.z;

				if(pm->x2 == -1)
				{
					// TODO: be a wall

				} else if(pm->x1 == cx2 && pm->z1 == cz2) {
					// 1 -> 2
					rx += (cx2-cx1) + (pm->x2-pm->x1);
					rz += (cz2-cz1) + (pm->z2-pm->z1);
					rcx = pm->x2;
					rcz = pm->z2;
					rot = (-pm->rot12) & 3;

				} else if(pm->x2 == cx2 && pm->z2 == cz2) {
					// 2 -> 1
					rx += (cx2-cx1) - (pm->x2-pm->x1);
					rz += (cz2-cz1) - (pm->z2-pm->z1);
					rcx = pm->x1;
					rcz = pm->z1;
					rot = pm->rot12 & 3;

				} else {
					// TODO: be a wall
				}

				float trx = rx;
				float trz = rz;
				float trvxx = rvxx;
				float trvxz = rvxz;
				float trvzx = rvzx;
				float trvzz = rvzz;

				switch(rot)
				{
					case 0:
						break;

					case 1:
						rx = (rcx + 0.5f) + (trz-(rcz + 0.5f));
						rz = (rcz + 0.5f) - (trx-(rcx + 0.5f));
						rvxx =  trvxz;
						rvxz = -trvxx;
						rvzx =  trvzz;
						rvzz = -trvzx;
						break;

					case 2:
						rx = (rcx + 0.5f)*2.0f - rx;
						rz = (rcz + 0.5f)*2.0f - rz;
						rvxx = -trvxx;
						rvxz = -trvxz;
						rvzx = -trvzx;
						rvzz = -trvzz;
						break;

					case 3:
						rx = (rcx + 0.5f) - (trz-(rcz + 0.5f));
						rz = (rcz + 0.5f) + (trx-(rcx + 0.5f));
						rvxx = -trvxz;
						rvxz =  trvxx;
						rvzx = -trvzz;
						rvzz =  trvzx;
						break;
				}

				cam.v.w.v.x = rx;
				cam.v.w.v.z = rz;
				cam.v.x.v.x = rvxx;
				cam.v.x.v.z = rvxz;
				cam.v.z.v.x = rvzx;
				cam.v.z.v.z = rvzz;
			}
		}
	}
}

int main(int argc, char *argv[])
{
	printf("Starting!\n");
	SDL_Init(SDL_INIT_VIDEO | SDL_INIT_NOPARACHUTE);

#ifndef WIN32
	signal(SIGINT,  SIG_DFL);
	signal(SIGTERM, SIG_DFL);
#endif

	SDL_WM_SetCaption("7DFPS 2014 (fanzyflani's entry)", NULL);
	screen = SDL_SetVideoMode(DEF_WIDTH, DEF_HEIGHT, 32, 0);
	sbuf = malloc(sizeof(uint32_t)*rwidth*rheight);
	tsbuf = malloc(sizeof(uint32_t)*rwidth*rheight);
	zbuf = malloc(sizeof(float)*rwidth*rheight);

	return mainloop();
}

