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

n, m = map(int, input().split())
imos = [0] * (n + 1)
for _ in range(m):
    l, r = map(int, input().split())
    l -= 1
    imos[l] += 1
    imos[r] -= 1
for i in range(1, n + 1):
    imos[i] += imos[i - 1]
ans = 1e9
for i in range(n):
    ans = min(ans, imos[i])
print(ans)











