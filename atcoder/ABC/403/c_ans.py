import io
import sys
from collections import defaultdict,deque,Counter
import math
# from more_itertools import chunked

# 下記に標準入力を記載
_INPUT = """\
5 5 10
2 2
3 4 4
1 1 1
1 4 1
1 4 2
1 4 4
1 2 4
3 3 2
3 5 4
3 2 1

"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

N, M, Q = map(int, input().split())
#setを使うことで計算量O(QM)をO(QlogM)にする
can_view = [set() for _ in range(N)]
#fxで状態管理し、計算量O(1)にする
can_view_all = [False] * N
for _ in range(Q):
    t, *q = map(int, input().split())
    x = q[0] - 1
    if t == 1:
        y = q[1]
        can_view[x].add(y)
        # print(can_view)
    elif t == 2:
        can_view_all[x] = True
    elif t == 3:
        y = q[1]
        print("Yes" if can_view_all[x] or y in can_view[x] else "No")
