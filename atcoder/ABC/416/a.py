import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np

# 下記に標準入力を記載
_InPUT = """\
10 6 8
xoxxooooxo


"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n, l, r = map(int, input().split())
s = input().strip()

cnt = 0

for index, s_item in enumerate(s):
    if index >= l - 1 and index <= r - 1:
        if s_item == 'o':
            cnt += 1
            if  cnt == r - l + 1:
                print('Yes')
        elif s_item == 'x':
            print('No')
            exit()
