static inline void trace_ray_prelude(level *lv, int cy, const vec4 *from, const vec4 *rayb, const vec4 *rdx, const vec4 *rdy, int dimx, int dimy)
{
	// work around a libgomp bug by shunting this out --GM
	int y, sy, cx, x;

	for(cx = 0; cx < dimx; cx += 32)
	for(sy = 0; sy < 32 && sy + cy < dimy; sy++)
	{
		y = sy + cy;

		vec4 rayl;
		rayl.m = _mm_add_ps(rayb->m, _mm_mul_ps(_mm_set1_ps(y), rdy->m));
		rayl.m = _mm_add_ps(rayl.m, _mm_mul_ps(_mm_set1_ps(cx), rdx->m));
		uint32_t *p = sbuf + y*rwidth + cx;
		float *dist = zbuf + y*rwidth + cx;
		for(x = 0; x < 32 && x + cx < dimx; x++)
		{
			rayl.m = _mm_add_ps(rayl.m, rdx->m);
			uint32_t seed = (x+cx) + y*y*(rwidth+1);
			seed *= seed * seed;
			seed *= seed * seed;
			*(p++) = col_ftoint(
				trace_ray(0, &seed, lv, dist++, from, &rayl,
				_mm_setr_ps(1.0f, 1.0f, 1.0f, 1.0f)));
		}
	}

}

// [p1,p2)
void trace_screen_centred(level *lv, int x1, int y1, int x2, int y2, const mat4 *cam)
{
	vec4 rayb, rdx, rdy;

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
#if !(defined(__GNUC__))
  #pragma omp parallel for
#endif
	for(cy = 0; cy < dimy; cy += 32)
	{
		trace_ray_prelude(lv, cy, from, &rayb, &rdx, &rdy, dimx, dimy);
	}

#ifdef POSTPROC_BLUR
	int blurctr;
	for(blurctr = 0; blurctr < POSTPROC_BLUR; blurctr++)
	{
		// Apply focal blur
		// TODO: isolate this to the camera region
		memcpy(tsbuf, sbuf, sizeof(uint32_t)*rwidth*rheight);

#if !(defined(__GNUC__))
  #pragma omp parallel for
#endif
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

