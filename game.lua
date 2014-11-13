oball = obj_new()
obx, oby, obz = 9.5, 0.3, 5.5
obj_set(oball, "sphere", 0.3, 0.9,
	obx, oby, obz,
	1.0, 1.0, 1.0)

obvx, obvz = 1, 0

function on_tick(sec_current, sec_delta)
	--
	obj_set(oball, "sphere", 0.3, 0.9,
		obx, oby, obz,
		1.0, 1.0, 1.0)
	
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

