import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import combinations,permutations

# 下記に標準入力を記載
_InPUT = """\
0 0 4 2
3 2 1
R 2
D 1
U 3
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# sys.setrecursionlimit(1000000)

r_t,c_t,r_a,c_a = map(int,input().split())
n,m,l = map(int,input().split())
for _ in range(m):
    s,a = input().split()
for _ in range(l):
    t,b = input().split()