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
for x in range(n, -1, -1):
    count = 0
    for a_i in a:
        if a_i >= x:
            count += 1
    if count >= x:
        print(x)
        break
