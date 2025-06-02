from collections import defaultdict
import io
import sys

# 下記に標準入力を記載
_INPUT = """\
7
3 1 4 1 5 9 2

"""


sys.stdin = io.StringIO(_INPUT)

A= int(input())
B= list(map(int, input().split()))
num = 0
for index,odd in enumerate(B):
    if (index+1) % 2 == 1:
        num += B[index]

print(num)

