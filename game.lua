local olist = {}
num = 10
for i=1,num do
	olist[i] = obj_new()
	obj_set(olist[i], "sphere", 0.5, 0.9,
		9.1 + 2.0*math.sin(i*math.pi*2/num),
		0.4,
		5.1 + 2.0*math.cos(i*math.pi*2/num),

		1.0, 1.0, 1.0)
end

