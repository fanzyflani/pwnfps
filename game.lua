oball = {}
opos = {
	{0.0, 0.0, 0.0, 0.3, 0.8, 0.8, 0.8, 0.6},

	{0.0, 0.3, 0.0, 0.1, 0.4, 0.4, 0.4, 0.2},

	{ 0.3, 0.0,  0.0, 0.1, 0.7, 0.7, 1.0, 0.4},
	{ 0.0, 0.0,  0.3, 0.1, 0.7, 1.0, 0.7, 0.4},
	{-0.3, 0.0,  0.0, 0.1, 1.0, 0.7, 0.7, 0.4},
	{ 0.0, 0.0, -0.3, 0.1, 0.5, 1.0, 1.0, 0.4},

	{ 0.3, 0.0,  0.1, 0.03, 0.4, 0.4, 0.4, 0.2},
	{ 0.1, 0.0,  0.3, 0.03, 0.4, 0.4, 0.4, 0.2},
	{-0.3, 0.0,  0.1, 0.03, 0.4, 0.4, 0.4, 0.2},
	{ 0.1, 0.0, -0.3, 0.03, 0.4, 0.4, 0.4, 0.2},
	{ 0.3, 0.0, -0.1, 0.03, 0.4, 0.4, 0.4, 0.2},
	{-0.1, 0.0,  0.3, 0.03, 0.4, 0.4, 0.4, 0.2},
	{-0.3, 0.0, -0.1, 0.03, 0.4, 0.4, 0.4, 0.2},
	{-0.1, 0.0, -0.3, 0.03, 0.4, 0.4, 0.4, 0.2},
}

obx, oby, obz = 9.5, 0.3, 5.5
obvx, obvz = 1, 0

for i=1,#opos do
	oball[i] = obj_new()
	obj_set(oball[i], "sphere", opos[i][4], opos[i][8],
		obx + opos[i][1], oby + opos[i][2], obz + opos[i][3],
		opos[i][5], opos[i][6], opos[i][7])
end


function on_tick(sec_current, sec_delta)
	--

	if math.fmod(sec_current, 0.5) < 0.15 then
		opos[2][5], opos[2][6], opos[2][7] = 0.3, 0.3, 1.3
	else
		opos[2][5], opos[2][6], opos[2][7] = 0.3, 0.3, 0.3
	end

	for i=1,#opos do
		local rx, ry, rz
		rx = opos[i][1]
		ry = opos[i][2]
		rz = opos[i][3]

		local rs = math.sin(sec_current*math.pi*2/2)
		local rc = math.cos(sec_current*math.pi*2/2)

		rx, rz = rc*rx + rs*rz, rc*rz - rs*rx

		obj_set(oball[i], "sphere", opos[i][4], opos[i][8],
			obx + rx, 
			oby + ry, 
			obz + rz, 
			opos[i][5], opos[i][6], opos[i][7])
	end
	
	local c1 = level_get(obx, obz)

	spd = 2.0
	nobx = obx + obvx*sec_delta*spd
	nobz = obz + obvz*sec_delta*spd

	local c2 = level_get(nobx + obvx*0.5, nobz + obvz*0.5)

	if c1 ~= c2 then
		if c2 == "." or ((c1 == "#" or c1 == "&") and c2 == '"') then
			obvx, obvz = obvz, -obvx
			nobx = obx + obvx*sec_delta*spd
			nobz = obz + obvz*sec_delta*spd
			c2 = level_get(nobx + obvx*0.5, nobz + obvz*0.5)
			if c2 == "." or ((c1 == "#" or c1 == "&") and c2 == '"') then
				obvx, obvz = -obvx, -obvz
				nobx = obx + obvx*sec_delta*spd
				nobz = obz + obvz*sec_delta*spd
				c2 = level_get(nobx + obvx*0.5, nobz + obvz*0.5)
			end
		end
	end
	obx, obz = nobx, nobz
	

	--print(sec_current, sec_delta)
end

