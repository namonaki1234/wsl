import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import islice

# 下記に標準入力を記載
_InPUT = """\
10 4
1 6
4 5
5 10
7 10

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n,m = map(int, input().split())
LR = [list(map(int, input().split())) for _ in range(m)]
# print(LR)

LR_counter = Counter(i for i in range(1, n + 1))
# print(LR_counter)

for l, r in LR:
    for i in range(l, r + 1):
        LR_counter[i] += 1
# LR_counter = sorted(LR_counter.items(), key=lambda x: x[1])
print(LR_counter)
LR_min = min(LR_counter.values())
print(LR_min-1)

# print(LR_counter[0][1]-1)