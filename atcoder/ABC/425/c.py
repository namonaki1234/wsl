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

n, q = map(int,input().split())
# a = deque(map(int,input().split()))
a = list(map(int,input().split()))
queries = [list(map(int,input().split())) for _ in range(q)]

for query in queries:
    if query[0] == 1:
        c = query[1]
        # a.rotate(-c)
        a = a[c:] + a[:c]

        # print(*a)
    elif query[0] == 2:
        l,r = query[1], query[2]
        # a = list(a)
        ans = 0
        for i in range(l,r+1):
            ans += a[i-1]
        # a = deque(a)
        print(ans)
