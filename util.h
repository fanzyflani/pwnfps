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

