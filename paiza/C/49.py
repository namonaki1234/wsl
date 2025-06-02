import io
import sys
import numpy as np

# 下記に標準入力を記載
_INPUT = """\
8
17
28
11
62
64
4
7
17
"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載 出力7

N = int(input())
f = [int(input()) for _ in range(N)]


count = 0
floor = 1
diff= f[0] - floor
for i in range(N-1):
    diff += abs(f[i+1] -f[i] )


print(diff)














