import io
import sys
from collections import defaultdict,deque,Counter
import math
# from more_itertools import chunked

# 下記に標準入力を記載
_INPUT = """\
1000000 500000
"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

n, k = map(int, input().split())
s = k
a = [1 for i in range(n + 1)]
for i in range(k, n + 1):
    a[i] = s
    s -= a[i-k]
    s += a[i]
    s %= 1000000000
print(a[n])
