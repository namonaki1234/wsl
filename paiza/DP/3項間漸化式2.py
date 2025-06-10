import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np

# 下記に標準入力を記載
_InPUT = """\
5
1
2
3
4
3
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

q = int(input())
k = [0] * q
for i in range(q):
    k[i] = int(input())

n = max(k)
a = [0] * (n + 1)
a[0] = 1
a[1] = 1
for i in range(2, n + 1):
    a[i] = a[i - 1] + a[i - 2]

for i in range(q):
    print(a[k[i]-1])