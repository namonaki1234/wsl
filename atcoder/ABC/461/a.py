import io
import sys

# 下記に標準入力を記載
_InPUT = """\
5 5
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

a,d = map(int,input().split())

if a <= d :
    print("Yes")
else:
    print("No")
