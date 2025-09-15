import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import combinations,permutations

# 下記に標準入力を記載
_InPUT = """\
4 10
30 300 3
60 45 4
90 45 5
120 45 2
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# sys.setrecursionlimit(1000000)

n,k = map(int,input().split())
for i in range(n):
    a,b,c = map(int,input().split())


print(n,k,a,b,c)
