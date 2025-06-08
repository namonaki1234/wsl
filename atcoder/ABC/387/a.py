import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import combinations

# 下記に標準入力を記載
_InPUT = """\
20 25
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

a, b = map(int, input().split())

print(a**2 + b**2 + 2*a*b)
