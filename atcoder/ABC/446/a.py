import io
import sys

# 下記に標準入力を記載
_InPUT = """\
Fred
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

s = input()

ans = "Of"+s[0].lower()+s[1:]
print(ans)