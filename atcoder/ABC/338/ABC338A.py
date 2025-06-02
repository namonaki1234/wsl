import re

u = input()

s = []
i = 0
while i < len(u) :
    s.append(u[i])
    i = i +1

m = re.findall(r'[A-Z]',u)

if len(m) == 1 :
    if m[0] == s[0]:
        print('Yes')
    else :print('No')

else :print('No')