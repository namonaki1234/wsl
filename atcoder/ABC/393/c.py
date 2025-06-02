import io
import sys
import numpy as np
from collections import deque

# 下記に標準入力を記載
_INPUT = """\
6 10
6 2
4 1
5 1
6 6
5 3
5 1
1 4
6 4
4 2
5 6

"""
2
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

N,M = list(map(int, input().split()))
if M != 0:
    uv = [map(int, input().split()) for _ in range(M)]
    u,v = [list(i) for i in zip(*uv)]
    uv =[list(i) for i in uv]
    print(uv)
roop_count = 0
double_count = 0

for i in range(M):
    if u[i] == v[i]:
        roop_count += 1
    for j in range(i,M):
        if i != j and i < j and u[i] == v[j] and u[j] == v[i]:
            double_count += 1
        
        elif i != j and i < j and u[i] == u[j] and v[j] == v[i]:
            double_count += 1
            
print(roop_count+double_count)



        


















