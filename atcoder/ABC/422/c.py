import io
import sys
from collections import defaultdict,deque,Counter

# 下記に標準入力を記載
_InPUT = """\
5
3 1 2
100 0 0
1000000 1000000 1000000
31 41 59
1000000000 10000 1
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# sys.setrecursionlimit(1000000)

t = int(input())
for _ in range(t):
    n_a,n_b,n_c = map(int,input().split())
    ans = min(n_a,n_b,n_c)

    print(ans)