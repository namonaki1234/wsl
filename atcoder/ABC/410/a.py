import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np

# 下記に標準入力を記載
_InPUT = """\
15
18 89 31 2 15 93 64 78 58 19 79 59 24 50 30
38

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())
a = list(map(int, input().split()))
k = int(input())

a.sort(reverse=True)

cnt = 0
for i in range(n):
    if a[i] >= k:
        cnt += 1

print(cnt)