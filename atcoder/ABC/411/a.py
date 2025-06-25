import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np

# 下記に標準入力を記載
_InPUT = """\
atcoder
7

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

p = input()
l = int(input())

if len(p) >= l:
    print("Yes")
else:
    print("No")
    exit()
