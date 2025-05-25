import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import islice

# 下記に標準入力を記載
_InPUT = """\
2025524202552420255242025524

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

s = input()

# Aの回数をカウント
count = len(s)

diff =[]
for i in range(len(s)):
    if i+1 < len(s):
        if s[i] < s[i+1]:
            diff.append(int(s[i])+10- int(s[i+1]))
        else:
            diff.append(int(s[i]) - int(s[i+1]))

print(sum(diff)+count+int(s[len(s)-1]))