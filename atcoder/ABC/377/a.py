import io
import sys

# 下記に標準入力を記載
_InPUT = """\
ARC
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

s = input()
s = list(s)
s.sort()

if s == ["A", "B", "C"]:
    print("Yes")
else:
    print("No")