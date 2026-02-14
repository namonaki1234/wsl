import io
import sys

# 下記に標準入力を記載
_InPUT = """\
rule
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

s = input()

if s[0] == s[len(s)-1]:
    print("Yes")
else:
    print("No")

