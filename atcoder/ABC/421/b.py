import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import combinations

# 下記に標準入力を記載
_InPUT = """\
90701 90204


"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

x,y = map(int, input().split())

next_a = [x,y]

def rev(x,y):
    z = x + y
    z = str(z)
    z = reversed(z)
    return int("".join(z))

for i in range(8):
    next_a.append(rev(next_a[i], next_a[i+1]))

# print(next_a)
print(next_a[-1])