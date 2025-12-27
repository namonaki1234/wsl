import io
import sys
from collections import Counter

# 下記に標準入力を記載
_InPUT = """\
8 3
20251227
438
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# 長さnの数字列s,長さmの数字列t

n,m = map(int,input().split())

s = input()
t = input()

min_diff = 100000
for i in range(n-(m-1)):
    target = s[i:i+m]
    # print(target)
    diff_sum = 0
    for j in range(m):
        diff_item = ((int(target[j])-int(t[j]))%10)
        # diff_item = ((int(target[j])-int(t[j]))%10)
        diff_sum += diff_item
        # diff_sum += diff_item if (10-diff_item) >= diff_item else 10-diff_item
        # print(diff_item)
    # diff_sum = ((int(target[0])-int(t[0]))%10)+((int(target[1])-int(t[1]))%10)
    # print(diff_sum)
    min_diff = min(diff_sum,min_diff)

print(min_diff)