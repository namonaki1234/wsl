import io
import sys

# 下記に標準入力を記載
_InPUT = """\
14 79
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

l,r = input().split()
# l,r = map(int(),input().split())
# print(l,r)
print(int(r)-int(l)+1)
