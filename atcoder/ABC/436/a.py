import io
import sys


# 下記に標準入力を記載
_InPUT = """\
12
vgxgpuam
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())
s = input()
ans = list(s)

for i in range(n-len(s)):
    ans.insert(0,"o")

print(''.join(ans))
