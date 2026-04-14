import io
import sys

# 下記に標準入力を記載
_InPUT = """\
5
ooooo
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())
s = input()
if "o" == s[0]:
    s = s.lstrip("o")

print(s)
