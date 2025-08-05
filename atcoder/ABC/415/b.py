import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
import itertools as it

# 下記に標準入力を記載
_InPUT = """\
.#.##..##.#.###....#


"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

s = input().strip()

cnt = 0
index = []
for i,item in enumerate(s):
    if item == "#":
        cnt += 1
        index.append(i+1)

    if cnt == 2:
        print(",".join(map(str,index)))
        cnt = 0
        index = []
