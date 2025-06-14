import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import combinations

# 下記に標準入力を記載
_InPUT = """\
5 6
4 3 1 2
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n, l = map(int, input().split())
if l % 3 != 0:
    print(0)
    exit()
d = list(map(int, input().split()))
x = 0
cnt = [0] * l
# xはmod6で累積和を表している
for i in range(n):
    if i != 0:
        x += d[i - 1]
    # xをmod6にしたうえで、位置のカウントを行っている
    x %= l
    cnt[x] += 1

ans = 0
for i in range(l // 3):
    ans += cnt[i] * cnt[i + l // 3] * cnt[i + 2 * l // 3]
print(ans)
