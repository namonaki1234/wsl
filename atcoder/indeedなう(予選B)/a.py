import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import combinations

# 下記に標準入力を記載
_InPUT = """\
3 3
2 5

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

x1, y1 = map(int, input().split())
x2, y2 = map(int, input().split())

intersection = abs(x2 - x1) + abs(y2 - y1) + 1

print(intersection)
