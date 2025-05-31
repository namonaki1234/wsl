import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import islice

# 下記に標準入力を記載
_InPUT = """\
5 2
1 2
3 4


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
        if LR_counter[i] == 1:
          del  LR_counter[i]
        if LR_counter[i] > 2:
            del LR_counter[i]
        # print(LR_counter)
# print(LR_counter)
LR_min = min(LR_counter.values())
print(LR_min-1)

# print(LR_counter[0][1]-1)