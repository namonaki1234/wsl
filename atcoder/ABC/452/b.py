import io
import sys
from collections import Counter

# 下記に標準入力を記載
_InPUT = """\
5 6
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

h,w = map(int, input().split())

start_end = "#" * w
middle = "#" + "." * (w - 2) + "#"
print(start_end)
for i in range(h - 2):
    print(middle)
print(start_end)



