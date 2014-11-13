int lf_obj_new(lua_State *L)
{
	part *pt = level_obj_new(lvroot);
	if(pt == NULL) return luaL_error(L, "obj_new: could not allocate object");

	lua_pushlightuserdata(L, pt);
	return 1;
}

int lf_obj_set(lua_State *L)
{
	part *pt = lua_touserdata(L, 1);
	if(pt == NULL) return luaL_error(L, "obj_set: pt cannot be nil");

	const char *typ = lua_tostring(L, 2);
	if(typ == NULL) return luaL_error(L, "obj_set: typ cannot be nil");

	if(!strcasecmp(typ, "sphere"))
	{
		pt->typ = P_SPHERE;

		pt->sph.r = lua_tonumber(L, 3);
		pt->sph.refl = lua_tonumber(L, 4);

		pt->sph.pos.v.x = lua_tonumber(L, 5);
		pt->sph.pos.v.y = lua_tonumber(L, 6);
		pt->sph.pos.v.z = lua_tonumber(L, 7);
		pt->sph.pos.v.w = 1.0f;

		pt->sph.col.c.b = lua_tonumber(L, 8);
		pt->sph.col.c.g = lua_tonumber(L, 9);
		pt->sph.col.c.r = lua_tonumber(L, 10);

	} else {
		return luaL_error(L, "obj_set: invalid typ \"%s\"", typ);
	}

	lua_pushlightuserdata(L, pt);
	return 1;
}

int lf_obj_free(lua_State *L)
{
	// TODO!
	part *pt = lua_touserdata(L, 1);
	if(pt == NULL) return luaL_error(L, "obj_set: pt cannot be nil");

	pt->typ = P_FREE;

	return 0;
}

int lf_level_get(lua_State *L)
{
	// TODO!
	return 0;
}

int lf_level_set(lua_State *L)
{
	// TODO!
	return 0;
}

lua_State *script_newvm(void)
{
	lua_State *L = luaL_newstate();

	// WARNING: be wary if someone gives you a random game.lua file
	// - check to see if it does any dodgy "io" calls and whatnot!
	luaL_openlibs(L);

	lua_pushcfunction(L, lf_obj_new); lua_setglobal(L, "obj_new");
	lua_pushcfunction(L, lf_obj_set); lua_setglobal(L, "obj_set");
	lua_pushcfunction(L, lf_obj_free); lua_setglobal(L, "obj_free");
	lua_pushcfunction(L, lf_level_get); lua_getglobal(L, "level_get");
	lua_pushcfunction(L, lf_level_set); lua_setglobal(L, "level_set");

	printf("Loading game.lua\n");
	if(luaL_loadfile(L, "game.lua") != 0)
	{
		printf("ERROR: game.lua failed to load (%s)\n", lua_tostring(L, -1));
		return NULL;
	}

	printf("Running game.lua\n");
	if(lua_pcall(L, 0, 0, 0) != 0)
	{
		printf("ERROR: game.lua failed to run (%s)\n", lua_tostring(L, -1));
		return NULL;
	}

	printf("Executed! Let's go!\n");

	return L;
}

