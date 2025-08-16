import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import combinations,permutations

# 下記に標準入力を記載
_InPUT = """\
5 3
13 13 13 13 2
5
12
13

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# sys.setrecursionlimit(1000000)

n,q = map(int, input().split())
a = list(map(int, input().split()))
b = [int(input()) for _ in range(q)]

# print(a)
# print(b)

for j in range(q):
    if b[j] > max(a):
        print(-1)
    elif b[j] == 1:
        print(1)
    elif b[j] == max(a):
        print(sum(a))
    elif 1 < b[j] < max(a):
        a_minus = [n-b[j]+1 for n in a if n-b[j]+1 < 0]
        print(sum(a_minus)+(b[j]-1)*n+1)

