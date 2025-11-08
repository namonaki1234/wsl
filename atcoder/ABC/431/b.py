import io
import sys

# 下記に標準入力を記載
_InPUT = """\
31
4
15 92 65 35
4
3
1
4
1
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# pは種類の番号,nは重りの数,xは最初の重さ

x = int(input())
n = int(input())
w = list(map(int, input().split()))
q = int(input())
cnt = [0] * n
for _ in range(q):
    p = int(input())
    if cnt[p - 1] >= 1:
        x -= w[p - 1]
        if x < 0:
            x = 0
        cnt[p - 1] -= 1
    else:
        x += w[p - 1]
        cnt[p - 1] += 1
    print(x)