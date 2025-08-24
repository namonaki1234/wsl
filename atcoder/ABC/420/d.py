import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import combinations,permutations

# 下記に標準入力を記載
_InPUT = """\
2 4
S.xG
#?o.
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# sys.setrecursionlimit(1000000)

h,w = map(int, input().split())
a = [list(input().strip()) for _ in range(h)]

s_pos = g_pos = None

for i, row in enumerate(a):
    for j, ch in enumerate(row):
        c = ch.lower()
        if c == 's':
            s_pos = (i, j)
        elif c == 'g':
            g_pos = (i, j)


def explore(start_row, start_col):
    maze_count = [[-1] * w for _ in range(h)]
    maze_count[start_row][start_col] = 0

    que = deque()
    que.append([start_row, start_col])

    while 0 < len(que):
        row, col = que.popleft()
        for dr, dc in [[1,0],[-1,0],[0,1],[0,-1]]:
            new_row = row + dr
            new_col = col + dc
            if 0<=new_row<h and 0<=new_col<w:
                if a[new_row][new_col] != '#' and maze_count[new_row][new_col] == -1:
                    maze_count[new_row][new_col] = maze_count[row][col] + 1
                    que.append([new_row, new_col])

    return maze_count

ans = explore(s_pos[0], s_pos[1])

print(ans)