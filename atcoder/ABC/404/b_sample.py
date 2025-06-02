import io
import sys
from collections import defaultdict,deque,Counter

# 下記に標準入力を記載
_INPUT = """\
3 5
.....
.#.#.
.....

"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

H, W = map(int, input().split())
B = [list(input()) for _ in range(H)]

for y in range(H):
    for x in range(W):
        if B[y][x] == '.':
            c = 0
            for dy in [-1, 0, 1]:
                for dx in [-1, 0, 1]:
                    if dx == 0 and dy == 0:
                        continue
                    ny = y + dy
                    nx = x + dx
                    if 0 <= ny < H and 0 <= nx < W:
                        if B[ny][nx] == '#':
                            c += 1
            B[y][x] = str(c)

for row in B:
    print(''.join(row))



    


    











