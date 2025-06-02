import io
import sys
# from more_itertools import chunked

# 下記に標準入力を記載
_INPUT = """\
5
3 1 4 1 5

"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

n = int(input())
a = list(map(int, input().split()))

suml = [0] * n
sumr = [0] * n
vis = [0] * (n + 1)

now = 0
for i in range(n):
    if vis[a[i]] == 0:
        now += 1
    vis[a[i]] = 1
    suml[i] = now

now = 0
vis = [0] * (n + 1)
for i in range(n - 1, -1, -1):
    if vis[a[i]] == 0:
        now += 1
    vis[a[i]] = 1
    sumr[i] = now

ans = 0
for i in range(n - 1):
    ans = max(ans, suml[i] + sumr[i + 1])

print(ans)
