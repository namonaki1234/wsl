import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import islice

# 下記に標準入力を記載
_InPUT = """\
12
11 3 8 9 5 2 10 4 1 6 12 7


"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())
p = list(map(int, input().split())) 

x =[]
y =[]

for i in range(n-1):
    if 2 <= i:
        if p[i-1]< p[i] > p[i+1]:
            x.append(i)
        elif p[i-1]> p[i] < p[i+1]:
            y.append(i)

ans = 0

if len(x) > 2 and len(y) > 2:
    for i in range(1,n-2):
        if y[i-1] < x[i] < y[i+1] and x[i-1] < y[i] < x[i+1]:
            if x[i] < y[i]:
                for j in range(x[i]):
                    if p[x[i]-j-1] < p[x[i]-j]:
                        ans += 1
            elif y[i] < x[i]:
                for j in range(y[i]):
                    if p[y[i]-j-1] < p[y[i]-j]:
                        ans += 1
                        
if (0 < len(x) < 3 and 0 < len(y) < 3) or (len(x) == 2 and 0 < len(y) <= 3) or (0 < len(x) <= 3 and len(y) == 2):
    for j in range(len(x)):
        if j!= 2 and (x[j] < y[j]):
            for k in range(x[j]):
                if p[x[j]-k-1] < p[x[j]-k]:
                    ans += 1
    
    for m in range(len(y)):
        if m!= 3 and (x[j] > y[j]):
            for l in range(y[j]):
                if p[y[j]-l-1] < p[y[j]-l]:
                    ans += 1

print(ans)
    











