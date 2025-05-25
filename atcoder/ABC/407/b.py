import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np

# 下記に標準入力を記載
_InPUT = """\
10 3


"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

x,y = map(int, input().split())

a,b =[i for i in range(1, 7) ], [i for i in range(1, 7) ]

complementary_event_count = 0

for i in a:
    for j in b:
        if (i+j) < x and abs(i-j) < y:
            complementary_event_count += 1

complementary_event_probability = complementary_event_count / (len(a) * len(b))
# print(complementary_event_count)
print(1-complementary_event_probability)
