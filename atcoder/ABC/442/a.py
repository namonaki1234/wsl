import io
import sys

# 下記に標準入力を記載
_InPUT = """\
jjjjjj
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

s = input()

count = 0
for s_i in s :
    if s_i == "i" or s_i == "j":
        count += 1

print(count)