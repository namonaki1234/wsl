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

memo_index = []
for i in range(len(s)):
    if i+1 == len(s):
        break
    if s[i] == s[i+1]:
        memo_index.append(i)


# print(memo_index)

ans = 0
if not memo_index:
    ans += len(s)*(len(s)+1) // 2
    print(ans%998244353)
    exit()

diff_index = [memo_index[0]]
for i in range(len(memo_index)-1):
    diff_index.append(memo_index[i+1]-memo_index[i])

# print(diff_index)

for i in range(len(diff_index)):
    if i == 0:
        ans += (diff_index[i]+2)*(diff_index[i]+2-1) // 2
        continue
    ans += (diff_index[i]+1)*(diff_index[i]+1-1) // 2

resume = len(s) - (memo_index[-1]+1)
ans += resume*(resume+1) // 2
print(int(ans)%998244353)