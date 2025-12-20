import io
import sys
from collections import Counter

# 下記に標準入力を記載
_InPUT = """\
3 4 5
12 3 5 7
6 10 11 9
1 2 4 8
2
4
9
6
11
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

h,w,n = map(int,input().split())

for i in range(h):
    a = list(map(int, input().split()))
