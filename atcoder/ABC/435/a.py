import io
import sys


# 下記に標準入力を記載
_InPUT = """\
29
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())

ans = int(0.5*(1+n)*n)

print(ans)
