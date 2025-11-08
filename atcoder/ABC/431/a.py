import io
import sys


# 下記に標準入力を記載
_InPUT = """\
1 1
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
h,b = map(int, input().split())

if h > b:
    print(h - b)
else:
    print(0)
