import io
import sys
import numpy as np
import pprint

# 下記に標準入力を記載
_INPUT = """\
5

"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載


N = int(input())
target = [["?"] * N for _ in range(N)]
# 格子を扱うときはなるべくi,jを0からN-1までの範囲で扱うと楽
for i in range(N):
    j = N - i - 1
    if i <= j:
        for x in range(i, j + 1):
            for y in range(i, j + 1):
                if i % 2 == 0:
                    target[x][y] = "#"
                else:
                    target[x][y] = "."
for row in target:
    print("".join(row))