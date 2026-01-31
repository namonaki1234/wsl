import io
import sys

# 下記に標準入力を記載
_InPUT = """\
beginner
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

s = input()

ans = s + "s"

print(ans)