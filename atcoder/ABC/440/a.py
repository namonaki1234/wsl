import io
import sys

# 下記に標準入力を記載
_InPUT = """\
432 1
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

x,y = map(int,input().split())

ans = x * (2 ** y)

print(ans)
