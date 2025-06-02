a,b = map(int,input().split())

c = []
for i in range(10):
    c.append(i)

d = a + b


c.remove(d)

print(c[0])
