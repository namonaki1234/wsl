import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
import itertools as it

# 下記に標準入力を記載
_InPUT = """\
6
g 4
j 1
m 4
e 4
d 3
i 4

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())
s = ""

for _ in range(n):
    c,l = input().split()
    l = int(l)
    if l > 100:
        print("Too Long")
        exit()
    s += c * l
    if len(s) > 100:
        print("Too Long")
        exit()

print(s)
