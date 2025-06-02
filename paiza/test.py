import io
import sys

# 下記に標準入力を記載
_INPUT = """\
1 2 5 3

"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載
a,b,c,d = map(int,input().split())
print((a+b)*(c-d))
print("Takahashi")