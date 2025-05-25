import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import islice

# 下記に標準入力を記載
_InPUT = """\
3 4
1 2 3 8
4 0 7 10
5 2 4 2
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

# XORスワップアルゴリズム
a = 10
b = 7
c = 0
d = 2
e = a^b^c^d
print(e)