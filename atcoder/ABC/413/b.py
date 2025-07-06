import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
import itertools as it

# 下記に標準入力を記載
_InPUT = """\
10
armiearggc
ukupaunpiy
cogzmjmiob
rtwbvmtruq
qapfzsitbl
vhkihnipny
ybonzypnsn
esxvgoudra
usngxmaqpt
yfseonwhgp

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())
s = [input().strip() for _ in range(n)]

sum_s = []
for i in range(n):
    for j in range(i + 1, n):
        sum_s.append(s[i] + s[j])
        sum_s.append(s[j] + s[i])
set_s = set(sum_s)
# print(set_s)
print(len(set_s))
