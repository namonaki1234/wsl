import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import combinations,permutations
from math import ceil

# 下記に標準入力を記載
_InPUT = """\
5 3
100 100 100 100 100
100 100 100 100 100
A 4 21
A 2 99
B 4 57



"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# sys.setrecursionlimit(1000000)

n,q = map(int, input().split())
a = list(map(int, input().split()))
b = list(map(int, input().split()))
queries = [input().split() for _ in range(q)]

for _ in range(q):
    min_sum_init = 0
    for k in range(n):
        min_sum_init += min(a[k], b[k])
# print(min_sum_init)

min_sum = min_sum_init
for query in queries:
    c,x,v = query
    x = int(x)
    v = int(v)
    if c == 'A':
        # min_sum = 0
        min_sum += (min(v, b[x-1]) - min(a[x-1], b[x-1]))
        # print(min(a[x-1], b[x-1]), min(v, b[x-1]))
        a[x-1] = v
        # for k in range(n):
            # print(v,b[k],a[k])
            # min_sum += min(a[k],b[k])
        print(min_sum)
    elif c == 'B':
        # min_sum = 0
        min_sum += (min(a[x-1], v) - min(a[x-1], b[x-1]))
        # print(min(a[x-1], b[x-1]), min(v, b[x-1]))
        b[x-1] = v
        # for k in range(n):
            # min_sum += min(a[k],b[k])
        print(min_sum)
