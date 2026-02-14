import io
import sys
from collections import Counter

# 下記に標準入力を記載
_InPUT = """\
6
abc
d
efghi
jkl
mnopq
r
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# 1<=n<=99
# longest = m 例1では11
# 奇数だから必ず先頭と末尾に同数の.があるはず
n = int(input())
# s =[0]*n
# len_s = [0]*n
arr = [{"s": 0, "len_s": 0} for _ in range(n)]
# for i in range(n):
#     s[i] = input()
#     len_s[i] = len(s[i])
for i in range(n):
    arr[i]["s"] = input()
    arr[i]["len_s"] = len(arr[i]["s"])

# print(s)
# print(len(s[-1]))

i_max = 0
max_len_s = arr[0]["len_s"]

for i in range(1, len(arr)):
    if arr[i]["len_s"] > max_len_s:
        max_len_s = arr[i]["len_s"]
        i_max = i

longest = arr[i_max]["s"]
# s_at_max = arr[i_max]["s"]

# print(i_max, max_len_s, s_at_max)
# longest = max()
ans = [0]*n
for j in range(n):
    add_s = "."*int((max_len_s-arr[j]["len_s"])/2)
    ans[j] = add_s + arr[j]["s"] + add_s

for k in range(n):
    print(ans[k])