import io
import sys

# 下記に標準入力を記載
_InPUT = """\
31415
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

s = list(input())
a_z_list = [chr(i) for i in range(ord('a'), ord('z') + 1)]

ans = ""
for s_i in s:
    if s_i not in a_z_list:
        ans += s_i

print(ans)
