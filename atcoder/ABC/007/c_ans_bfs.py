import io
import sys
from collections import deque
import itertools
from atcoder.dsu import DSU

# 再帰の上限を増やす（今回はBFSなので必須ではないが、そのままでもOK）
sys.setrecursionlimit(10**6)

# 下記に標準入力を記載
_InPUT = """\
7 8
2 2
4 5
########
#......#
#.######
#..#...#
#..##..#
##.....#
########
"""
sys.stdin = io.StringIO(_InPUT)

# ここからコードを記載

# 入力
r, c = map(int, input().split())
sx, sy = map(int, input().split())
tx, ty = map(int, input().split())

# 0-index に直す
sx -= 1
sy -= 1
tx -= 1
ty -= 1

grid = [input().strip() for _ in range(r)]

# 4方向移動（dx,dyの組み合わせで(x+1, y) → 下、(x, y+1) → 右、(x-1, y) → 上、(x, y-1) → 左）となっている
dx = [1, 0, -1, 0]
dy = [0, 1, 0, -1]

# 距離配列
dist = [[-1] * c for _ in range(r)]

# BFS
q = deque()
q.append((sx, sy))
dist[sx][sy] = 0

while q:
    x, y = q.popleft()

    for i in range(4):
        nx = x + dx[i]
        ny = y + dy[i]

        # 範囲外
        if nx < 0 or nx >= r or ny < 0 or ny >= c:
            continue

        # 壁
        if grid[nx][ny] == "#":
            continue

        # 訪問済み
        if dist[nx][ny] != -1:
            continue

        dist[nx][ny] = dist[x][y] + 1
        q.append((nx, ny))

print(dist[tx][ty])