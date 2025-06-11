import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import islice

# 下記に標準入力を記載
_InPUT = """\
5 10
2 5
1 5
1 2
2 4
2 2
5 5
2 4
1 2
2 2
2 3

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n,m = map(int, input().split())
LR = [list(map(int, input().split())) for _ in range(m)]
# print(LR)

canon = [0] * (n + 2)
for l, r in LR:
  canon[l] += 1
  canon[r+1] -= 1

ans = 1e9
for i in range(1, n + 1):
  canon[i] += canon[i - 1]
  ans = min(ans, canon[i])
print(ans)

# print(canon)