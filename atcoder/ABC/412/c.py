import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import combinations

# 下記に標準入力を記載
_InPUT = """\
3
4
1 3 2 5
2
1 100
10
298077099 766294630 440423914 59187620 725560241 585990757 965580536 623321126 550925214 917827435
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

t = int(input())
cases = []
for _ in range(t):
    n = int(input())
    s = list(map(int, input().split()))
    cases.append((n, s))

ans = []
for n, s in cases:
    s.sort()
    cnt = 0
    for i in range(n):
        for j in range(i + 1, n):
            if s[i] <= 2 * s[j]:
                cnt += 1
    if cnt == n - 1:
        ans.append(len(s))

for i in range(len(ans)):
        print(ans[i])
