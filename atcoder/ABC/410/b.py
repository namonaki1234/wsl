import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np

# 下記に標準入力を記載
_InPUT = """\
6 20
4 6 0 3 4 2 6 5 2 3 0 3 2 5 0 3 5 0 2 0
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n, q = map(int, input().split())
x = list(map(int, input().split()))

box = [0] * n

# for i in range(q):
#     if x[i] >= 1:
#         # box[i] = x[i]
#         box[i] += 1
#         print(x[i])
#     elif x[i] == 0:
#         # print(box.index(min(box)))
#         box[box.index(min(box))] += 1
#         print(box.index(min(box)))
#     print(box)

ans = []
for i in range(q):
    if x[i] >= 1:
        # print(x[i])
        ans.append(x[i])
        box[x[i]-1] += 1
    elif x[i] == 0:
        min_index = box.index(min(box))
        # print(min_index + 1)
        ans.append(min_index + 1)
        box[min_index] += 1

print(*ans)