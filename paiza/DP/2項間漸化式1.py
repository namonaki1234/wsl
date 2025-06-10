import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np

# 下記に標準入力を記載
_InPUT = """\
0 7 9
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

x, d, k = map(int, input().split())

a =[0] * (k + 1)
a[0] = x
for i in range(2, k + 1):
    a[i] = a[i - 1] + d

# print(*a, sep="\n")
print(a[k])