import re

s = input()
e = re.search(r'=+',s)

if  e:
    if s[0] == '<' and s[1:len(s)-1] == e.group() and s[len(s)-1] == '>':
        print('Yes')
    elif s[1:len(s)-1] != e.group() :
        print('No')
    else :
        print('No')
else :
    print('No')

