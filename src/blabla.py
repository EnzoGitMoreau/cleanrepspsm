l = []
for a in range(1,60):
	for b in range(1,60):
		for c in range(1,60):
			if(a+b+c == 60):
				l.append((a,b,c))

print(len(l))
