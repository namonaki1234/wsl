import io
import sys
from collections import Counter, deque
import itertools
import bisect
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

n = 0
l = 0
r = 0
# iteration = len(s)
while l < len(s):
    while r < len(s):
        target_count = r - l + 1
        if target_count %2 == 0:
            r+=1
            continue
        if dict_s[(l)+((target_count+1)//2 - 1)] == "C":
            count+=1
        r+=1
    l+=1
    r = l

    
# for i in range(len(s)):
#     for j in range(i+1,len(s)+1):
#         target_count = j-i
#         if target_count % 2 == 0:
#             continue
#         if dict_s[(i)+((target_count+1)//2 - 1)] == "C":
#             count+=1
print(count)
# 5

# S 中の # 1 文字目から # 5 文字目までを抜き出した ABCCA
# S 中の # 2 文字目から # 4 文字目までを抜き出した BCC
# S 中の # 3 文字目から # 3 文字目までを抜き出した C
# S 中の # 3 文字目から # 5 文字目までを抜き出した CCA
# S 中の # 4 文字目から # 4 文字目までを抜き出した C
