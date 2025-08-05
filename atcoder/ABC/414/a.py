import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np

# 下記に標準入力を記載
_InPUT = """\
10 8 14
5 20
14 21
9 21
5 23
8 10
0 14
3 8
2 6
0 16
5 20

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n,l,r = map(int, input().split())
x,y = [], []
for _ in range(n):
    a,b = map(int, input().split())
    x.append(a)
    y.append(b)

cnt = 0
for i in range(n):
    if x[i] <= l and y[i] >= r:
        cnt += 1
    continue

print(cnt)
