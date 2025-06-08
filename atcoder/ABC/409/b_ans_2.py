import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np

# 下記に標準入力を記載
_InPUT = """\
7
1 6 2 10 2 3 2
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())
a = list(map(int, input().split()))
a.sort(reverse=True)

# xは1-indexdで考える
for x in range(n, 0, -1):
    if a[x - 1] >= x:
        print(x)
        break
else:
    print(0)
