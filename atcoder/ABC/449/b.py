import io
import sys
from collections import Counter

# 下記に標準入力を記載
_InPUT = """\
7 9 5
2 4
1 3
2 1
2 1
1 3
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

h,w,q = map(int,input().split())
for i in range(q):
    query = list(map(int,input().split()))
    # print(query)
    if query[0] == 1:
        r = query[-1]
        print(w*r)
        h -= r
    elif query[0] == 2:
        c = query[-1]
        print(h*c)
        w -= c

