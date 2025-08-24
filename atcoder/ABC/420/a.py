import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np

# 下記に標準入力を記載
_InPUT = """\
12 12


"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

x, y = map(int, input().split())

after = (x + y) % 12 if (x + y) % 12 != 0 else 12

print(after)