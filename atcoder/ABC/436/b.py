import io
import sys
from collections import Counter

# 下記に標準入力を記載
_InPUT = """\
5
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())

cols = n;rows = n
grid = [[0 for _ in range(cols)] for _ in range(rows)]

r,c = 0,0
for i in range(n):
    for j in range(n):
        if i == j == 0:
            grid[0][int((n-1)/2)]=1
            r = 0;c = int((n-1)/2);k = 1
            continue

        serch_index = [(r-1)%n,(c+1)%n]
        if not grid[serch_index[0]][serch_index[1]]:
            grid[serch_index[0]][serch_index[1]] = k + 1
            r = serch_index[0];c = serch_index[1];k = k + 1
        else:
            grid[(r+1)%n][c] = k + 1
            r = (r+1)%n;c = c;k = k + 1

# print(grid)
for g in grid:
    print(*g)