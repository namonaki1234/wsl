import io
import sys
from collections import defaultdict,deque,Counter

# 下記に標準入力を記載
_InPUT = """\
5 7
1 2 4 8 16
2 1 5
1 4
1 5
2 1 5
2 2 4
1 1
2 3 3

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n, q = map(int, input().split())
a = list(map(int, input().split()))
b = a + a
for i in range(2 * n - 1, 0, -1):
    b[i - 1] += b[i]
rui_c = 0
for _ in range(q):
    query = list(map(int, input().split()))
    if query[0] == 1:
        c = query[1]
        rui_c += c
        rui_c %= n
    else:
        l, r = query[1] - 1, query[2]
        print(b[l + rui_c] - b[r + rui_c])
