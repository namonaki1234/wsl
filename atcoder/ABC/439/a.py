import io
import sys

# 下記に標準入力を記載
_InPUT = """\
11
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())

ans = 2**(n) - 2 * n

print(ans)
