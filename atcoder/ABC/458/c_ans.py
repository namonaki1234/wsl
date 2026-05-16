import io
import sys
from collections import Counter, deque
import itertools
import bisect
import re
# from atcoder.dsu import DSU

# 再帰の上限を増やす（Pythonでは必須のおまじない）
sys.setrecursionlimit(10**6)

# 下記に標準入力を記載
_InPUT = """\
ABCCA
"""
sys.stdin = io.StringIO(_InPUT)

# ここからコードを記載

s = input()
# print(s[1:2])
count = 0
# print(s[(len(s)+1)//2 - 1])
# s_q = deque(s)
dict_s = dict(enumerate(s))

# max_c = max((len(m) for m in re.findall(r"c+", s)), default=0)
# 一個の時
# count += s.count("C")

for i in range(len(s)):
    if s[i] == "C":
        count += min(i+1 - 1 + 1,len(s) - (i+1) + 1)
print(count)
# 5

# S 中の # 1 文字目から # 5 文字目までを抜き出した ABCCA
# S 中の # 2 文字目から # 4 文字目までを抜き出した BCC
# S 中の # 3 文字目から # 3 文字目までを抜き出した C
# S 中の # 3 文字目から # 5 文字目までを抜き出した CCA
# S 中の # 4 文字目から # 4 文字目までを抜き出した C
