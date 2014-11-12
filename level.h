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

