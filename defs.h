#define EPSILON 0.0000000000001f

// yes you have enough RAM
#define OBJ_MAX 10000

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
} __attribute__((aligned(16))) vec4;

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
	P_FREE,

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
} __attribute__((aligned(16)));

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
	part objs[OBJ_MAX];
	ssize_t objs_num;

	int sx, sz;

	portal pmap[26];

	char data[64][64];

	part **parts[64][64];
	uint16_t parts_num[64][64];
	uint16_t parts_max[64][64];

} level;

