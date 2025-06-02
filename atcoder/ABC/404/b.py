import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np

# 下記に標準入力を記載
_INPUT = """\
13
.#..###..##..
#.#.#..#.#.#.
#.#.###..#...
###.#..#.#.#.
#.#.###..##..
.............
..#...#....#.
.##..#.#..##.
#.#..#.#.#.#.
####.#.#.####
..#..#.#...#.
..#...#....#.
.............
.............
.#....#...#..
.#...#.#..#..
####.#.#.####
.#.#.###..#.#
.##....#..##.
.#....#...#..
.............
..##..###.#.#
.#.#.#..#.###
.#.#..###.#.#
.#.#.#..#.#.#
..##..###..#.

"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

N = int(input())
S = [input() for _ in range(N)]
T = [input() for _ in range(N)]

t = np.array([list(row) for row in T])
s = np.array([list(row) for row in S])

def rotate(grid, times):
    return np.rot90(grid, -times)

min_ans = float('inf')
for rot in range(4):
    rot_s = rotate(s, rot)
    diff = (rot_s != t)
    changes = np.count_nonzero(diff)

    total = rot + changes
    min_ans = min(min_ans, total)

print(min_ans)

    


    











