import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np

# 下記に標準入力を記載
_InPUT = """\
20 4 8
abcdefghijklmnopqrst

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n, a, b = map(int, input().split())
s = list(input().strip())

ans = s[a:len(s) - b]
print(''.join(ans))
