import io
import sys
import numpy as np

# 下記に標準入力を記載
_INPUT = """\
4 3
2556 3424 77
137
721
984
999
"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

N,M = map(int,input().split())
A,B,C = map(int,input().split())
R = [int(input()) for _ in range(N)]

count = 0
for i in range(N):
    benefit = C*R[i]-A-B*M
    if benefit < 0:
        count += 1
    else:
        pass
print(count)












