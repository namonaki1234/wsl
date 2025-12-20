import io
import sys


# 下記に標準入力を記載
_InPUT = """\
8 0
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

a,b = map(int,input().split())

ans = a * 12 + b 

print(ans)
