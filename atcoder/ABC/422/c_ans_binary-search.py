import io
import sys

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

def solve():
    a, b, c = map(int, input().split())

    def f(m):
        if m > min(a, c):
            return False
        rem = (a - m) + (c - m) + b
        return rem >= m

    ng, ok = 10**10, -1
    while ng - ok > 1:
        mid = (ok + ng) // 2
        if f(mid):
            ok = mid
        else:
            ng = mid
    print(ok)


for _ in range(int(input())):
    solve()
