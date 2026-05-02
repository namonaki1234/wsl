import io
import sys
from collections import Counter, deque
import itertools
# from atcoder.dsu import DSU

# 再帰の上限を増やす（Pythonでは必須のおまじない）
sys.setrecursionlimit(10**6)

# 下記に標準入力を記載
_InPUT = """\
cabcabcbcaccacbcbcaabacbacaabccacbccbcacbacbacabcacabcaccaaaaabababcbabacaccabbcacbcbcbcababcbcbabca
"""
sys.stdin = io.StringIO(_InPUT)

# ここからコードを記載

s = input()

ord_a = ord("a")

# print(ord_a)

r = 0
l = 1

ans = len(s) # 文字一個ずつの分は最初に入れる

while (r < len(s) - 1) and (l < len(s) - 1):
    if s[r] == s[l]:
        continue
    if l == len(s) -1:
        r += 1
        l == 0
    ans += 1
    l += 1

# for i in range(len(s)):
    # for j in range(i+1,len(s)):
    #     if s[i] == s[j]:
    #         continue
    #     ans += 1

print(ans)