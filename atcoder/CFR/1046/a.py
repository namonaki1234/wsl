import io
import sys

# 下記に標準入力を記載
_InPUT = """\
11
1 4 1 4
4 1 4 1
1 4 2 5
0 100 0 100
1 4 2 9
3 1 13 5
8 11 17 36
19 41 30 50
20 38 30 60
0 0 0 0
100 100 100 100
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
def three_check(x, y):
    return x <= 2 * y + 2 and y <= 2 * x + 2

t = int(input())
for _ in range(t):
    a, b, c, d = map(int, input().split())
    # 前半: a, b
    # 後半: c-a, d-b が追加得点
    if three_check(a, b) and three_check(c - a, d - b):
        print("YES")
    else:
        print("NO")
