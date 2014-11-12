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

