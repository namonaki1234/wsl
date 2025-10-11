import io
import sys
from collections import defaultdict,deque,Counter

# 下記に標準入力を記載
_InPUT = """\
6
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())

a_1_5 = []
a = 1

# nが1~5の場合
for i in range(5):
    a_1_5.append(a)
    a += a

if 1 <= n <= 5:
        print(a_1_5[n-1])
        exit()

def f(n: int) -> int:
    return sum(int(digit) for digit in str(abs(n)))

# print(f(16))
a = 0
for i in range(6, n+1):
    a = (a_1_5[i-2]) + f(a_1_5[i-2])
    # print(a)
    a_1_5.append(a)


print(a)