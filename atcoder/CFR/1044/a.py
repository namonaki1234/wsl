import io
import sys
from collections import defaultdict,deque,Counter
# import numpy as np
from itertools import combinations,permutations

# 下記に標準入力を記載
_InPUT = """\
5
2
5 5
4
6 3 6 9
2
2 3
7
30 10 12 10 10 9 18
5
2 4 8 16 32

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# sys.setrecursionlimit(1000000)

t = int(input())
for _ in range(t):
    n = int(input())
    a = list(map(int, input().split()))

    set_a = set(a)

    if len(set_a) == len(a):
        print("No")
    else:
        print("Yes")
