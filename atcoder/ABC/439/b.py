import io
import sys
from collections import Counter

# 下記に標準入力を記載
_InPUT = """\
1111
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())

n_str = str(n)
iteration = len(n_str)

n_sum = 100
for i in range(100000):
    n_sum = 0
    for j in range(len(n_str)):
        n_sum += int(n_str[j])**2
    # print(n_sum)
    n_str = str(n_sum)
    if n_sum == 1 :
        print("Yes")
        exit()
    elif len(n_str) == 1:
        print("No")
        exit()
    n_sum = 0

print("No")