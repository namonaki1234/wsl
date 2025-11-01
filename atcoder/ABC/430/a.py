import io
import sys


# 下記に標準入力を記載
_InPUT = """\
10 20 30 40
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

a, b, c, d = map(int, input().split())

if c >= a:
    if d >= b:
        print("No")
    else:
        print("Yes")
else:
    print("No")
