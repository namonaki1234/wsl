import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np

# 下記に標準入力を記載
_InPUT = """\
1 10000
100

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n,m = map(int, input().split())
a = list(map(int, input().split()))

if np.sum(a) <= m:
    print("Yes")
else:
    print("No")