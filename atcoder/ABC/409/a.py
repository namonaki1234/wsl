import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np

# 下記に標準入力を記載
_InPUT = """\
10
xoooxoxxxo
ooxooooxoo

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())
t = input()
a = input()

for i in range(n):
    if t[i] == 'o' and a[i] == 'o':
        print('Yes')
        exit()

print('No')