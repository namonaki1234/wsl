import io
import sys

# 下記に標準入力を記載
_InPUT = """\
8 20
9 19 14 17 17 4 18 4
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# nは整数の長さ

n,x = map(int,input().split())
a = list(map(int,input().split()))
# print(a)

for i in range(n):
    if a[i] < x:
        x = a[i] 
        print(1)
    else:
        print(0)
