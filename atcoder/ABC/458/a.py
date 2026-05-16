import io
import sys

# 下記に標準入力を記載
_InPUT = """\
burger
1
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

s = input()
n = int(input())


print(s[n:-n])