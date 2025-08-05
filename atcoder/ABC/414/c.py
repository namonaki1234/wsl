import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import combinations

# 下記に標準入力を記載
_InPUT = """\
8
999999999999

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

a = int(input())
n = int(input())

ans = []
#n以下の値は10進法のとき回文となるか
for i in range(1,n + 1):
    if str(i) == str(i)[::-1]:
        #n以下の値はa進法のとき回文となるか
        s = ''
        x = i
        while x > 0:
            s += str(x % a)
            x //= a
        if s == s[::-1]:
            # print(i)  # 回文であれば出力
            ans.append(i)

print(sum(ans))  # 回文の合計を出力