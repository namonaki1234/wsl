import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import combinations,permutations
from math import ceil

# 下記に標準入力を記載
_InPUT = """\
6
91 999999986
53 999999997
32 999999932
14 999999909
49 999999985
28 999999926



"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# sys.setrecursionlimit(1000000)

n = int(input())
r = []
c = []

for _ in range(n):
    a, b = map(int, input().split())
    r.append(a)
    c.append(b)

# print(r)
# print(c)

r_min = min(r)
r_max = max(r)
c_min = min(c)
c_max = max(c)

r_median = (r_min + r_max) / 2
c_median = (c_min + c_max) / 2

# print(r_median, c_median)

r_diff = ceil(r_median - r_min)
c_diff = ceil(c_median - c_min)

print(max(r_diff, c_diff))
