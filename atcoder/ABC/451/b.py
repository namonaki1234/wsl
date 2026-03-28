import io
import sys
from collections import Counter

# 下記に標準入力を記載
_InPUT = """\
5 4
1 2
2 1
3 1
2 2
2 4
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n,m = map(int, input().split())

c = [list(map(int, input().split())) for _ in range(n)]

# print(c)

now = [0]*n
next = [0]*n
ans = [0]*n
for c_i in c:
    now[c_i[0]] += 1 
    next[c_i[-1]] += 1
# print(now,next)

for i in range(1,m+1):
    ans[i] = next[i] - now[i]
    print(ans[i])

