import io
import sys
from collections import Counter

# 下記に標準入力を記載
_InPUT = """\
5
8 6 10 5 7
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())
a = list(map(int, input().split()))

# sum_a = []
# sum_plus = 0
# for i in range(n):
#     sum_plus += a[i]
#     sum_a.append(sum_plus)

# print(sum_a)
