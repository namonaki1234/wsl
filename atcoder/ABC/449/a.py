import io
import sys

# 下記に標準入力を記載
_InPUT = """\
98
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

d = int(input())
pi = 3.141592653589793
ans = pi * (d/2)**2
print(ans)