import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import combinations

# 下記に標準入力を記載
_InPUT = """\
8 7
.######
##....#
#.###.#
#.#.#.#
#.#.#.#
#.#####
#...#..
#####..
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

h,w = map(int,input().split())
s = [input() for _ in range(h)]
start = (0,0)
goal = (h-1,w-1)
directions = [(1,0),(-1,0),(0,1),(0,-1)]
dist = np.full((h,w),-1)
dist[0][0] = 0
q = deque([start])
while q:
    v = q.popleft()
    if v == goal:
        break
    for d in directions:
        nv = (v[0] + d[0], v[1] + d[1])
        if 0 <= nv[0] < h and 0 <= nv[1] < w:
            if s[nv[0]][nv[1]] == "." and dist[nv[0]][nv[1]] == -1:
                dist[nv[0]][nv[1]] = dist[v[0]][v[1]] + 1
                q.append(nv)
if dist[goal[0]][goal[1]] == -1:
    print("No")
else:
    print(dist[goal[0]][goal[1]])