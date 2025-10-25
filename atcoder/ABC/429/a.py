import io
import sys


# 下記に標準入力を記載
_InPUT = """\
3 5
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n,m = map(int, input().split())

for i in range(1, n+1):
    if i <= m:
        print("OK")
    else:
        print("Too Many Requests")