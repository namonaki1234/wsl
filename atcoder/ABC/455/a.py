import io
import sys

# 下記に標準入力を記載
_InPUT = """\
1 3 7
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

a,b,c = input().split()

if a != b and b == c:
    print("Yes")
else:
    print("No")
