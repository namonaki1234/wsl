import io
import sys

# 下記に標準入力を記載
_InPUT = """\
460 33
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n,m = map(int,input().split())

cnt = 0
while m != 0:
    x = n % m
    m = x
    cnt += 1

print(cnt)