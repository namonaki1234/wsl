import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
import itertools as itbbbE

# 下記に標準入力を記載
_InPUT = """\
#..#.
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

# s = input().strip()
print(("#"+input()).replace("#.","#o")[1:])
