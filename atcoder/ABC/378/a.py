import io
import sys
from collections import Counter

# 下記に標準入力を記載
_InPUT = """\
1 1 1 1
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# aの値は色を示す

a_1, a_2, a_3, a_4 = map(int, input().split())

a = [a_1, a_2, a_3, a_4]

a_counter = Counter(a)
# print(a_counter.keys(), a_counter.values())

double_count = 0
for value in a_counter.values():
    if value == 4:
        double_count += 2
        continue
    if value < 2:
        continue
    else:
        double_count += 1
print(double_count)