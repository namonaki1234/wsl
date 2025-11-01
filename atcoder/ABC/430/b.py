import io
import sys

# 下記に標準入力を記載
_InPUT = """\
10 3
..#.......
.###......
.#.#......
#####.....
#...#.....
......####
......#..#
......#...
......#..#
......####
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# #が黒、.が白

n,m = map(int,input().split())
s_grid = [input() for _ in range(n)]
# print(s_grid)
mm_grid = []

for i in range(n-(m-1)):
    for j in range(n-(m-1)):
        mm = []
        for k in range(m):
            mm.append(s_grid[i+k][j:j+m])
        mm_grid.append(mm)
    
# print(mm_grid)
set_mm = set(tuple(x) for x in mm_grid)
# print(set_mm)
print(len(set_mm))