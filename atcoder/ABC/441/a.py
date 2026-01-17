import io
import sys

# 下記に標準入力を記載
_InPUT = """\
5 5
10 1000
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

p,q = map(int,input().split())
x,y = map(int,input().split())

if p <= x <= p+99 and q <= y <= q+99:
    print("Yes")
else:
    print("No")