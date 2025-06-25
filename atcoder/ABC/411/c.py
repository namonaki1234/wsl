import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import combinations

# 下記に標準入力を記載
_InPUT = """\
5 7
2 3 3 5 1 5 2
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n, q = map(int, input().split())
a = list(map(int, input().split()))

box = [0] * (n)
l = []
r = []
for q_item in q:
    box[q_item-1] += 1

