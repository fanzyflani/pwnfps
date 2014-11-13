#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <errno.h>

#include <assert.h>
#include <math.h>

#include <SDL.h>
#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>

#include <mmintrin.h>
#include <xmmintrin.h>
#include <emmintrin.h>

#ifndef WIN32
#include <signal.h>
#endif

//define inline 
#include "defs.h"

int rwidth = DEF_RWIDTH;
int rheight = DEF_RHEIGHT;
int rscale = DEF_SCALE;

SDL_Surface *screen = NULL;
uint32_t *sbuf = NULL;
uint32_t *tsbuf = NULL;
float *zbuf = NULL;
level *lvroot = NULL;

#include "util.h"
#include "trace.h"
#include "screen.h"
#include "level.h"
#include "script.h"

int mainloop(void)
{
	int i;

	SDL_Event ev;
	mat4 cam;
	vec4 gravity;

	// TODO: move this part to Lua
	lvroot = level_load("level.txt");
	level *lv = lvroot;

	lua_State *L = script_newvm();
	if(L == NULL)
		return 1;

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

	int sph_count = 0;
	part *sph[sph_count];
	for(i = 0; i < sph_count; i++)
	{
		sph[i] = level_obj_new(lv);
		sph[i]->typ = P_SPHERE;
		sph[i]->sph.r = 0.3f;
		sph[i]->sph.refl = 0.6f;
		sph[i]->sph.col.m = _mm_setr_ps(1.0f, 0.5f, 0.5f, 1.0f);
	}

	for(;;)
	{
		level_prepare_render(lv);

		for(i = 0; i < sph_count; i++)
		{
			sph[i]->sph.pos.m = _mm_setr_ps(
				2.0f+lv->sx + sinf(M_PI*(sec_current*0.2f + i*2/(float)sph_count))*3.7f,
				0.3f + 0.1f + sinf(2.0f*M_PI*(sec_current*0.2f + i*10/(float)sph_count))*0.1f,
				2.0f+lv->sz + cosf(M_PI*(sec_current*0.2f + i*2/(float)sph_count))*3.7f,
				1.0f);

			//level_part_add(lv, sph[i]);
		}
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

		lua_getglobal(L, "on_tick");
		if(lua_isnil(L, -1))
		{
			printf("ERROR: on_tick not defined\n");
			return 2;
		}

		lua_pushnumber(L, sec_current);
		lua_pushnumber(L, tdiff);
		if(lua_pcall(L, 2, 0, 0) != 0)
		{
			printf("ERROR: on_tick broke (%s)\n", lua_tostring(L, -1));
			return 2;
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
	sbuf = malloc(sizeof(uint32_t)*rwidth*rheight + 16);
	tsbuf = malloc(sizeof(uint32_t)*rwidth*rheight + 16);
	zbuf = malloc(sizeof(float)*rwidth*rheight + 16);
	sbuf += (((intptr_t)sbuf)&12)>>2;
	tsbuf += (((intptr_t)sbuf)&12)>>2;
	zbuf += (((intptr_t)sbuf)&12)>>2;

	return mainloop();
}

