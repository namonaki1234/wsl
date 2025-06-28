import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np

# 下記に標準入力を記載
_InPUT = """\
aBCdE
abcdcba


"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

s = input()
t = input()

res = True
n = len(s)

for i in range(1, n):
    if s[i].isupper():
        if s[i - 1] not in t:
            res = False

print("Yes" if res else "No")
