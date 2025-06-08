import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np

# 下記に標準入力を記載
_InPUT = """\
11

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

x = int(input())

sumproduct = []
for i in range(9):
    for j in range(9):
        sumproduct.append((i+1)*(j+1))

if x in sumproduct:
    Counter = Counter(sumproduct)
    # print(Counter[x])
    ans = 2025- Counter[x]*x
else:
    ans = 2025
print(ans)
