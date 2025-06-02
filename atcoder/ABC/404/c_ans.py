import io
import sys
from collections import defaultdict,deque,Counter
import math
# from more_itertools import chunked

# 下記に標準入力を記載
_INPUT = """\
5 4
2 1 2
3 3 4 5
3 1 2 5
1 3
1 3 2 5 4


"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

n, m = map(int, input().split())

a = []  # 各グループのリスト
idx = [[] for _ in range(n)]  # 各数字がどのグループに含まれるか
cnt = [0] * m

for i in range(m):
    parts = list(map(int, input().split()))
    k = parts[0]
    rest = parts[1:]
    if len(rest) < k:
        rest += list(map(int, input().split()))
    cnt[i] = k
    a.append([x - 1 for x in rest])  # 0-indexed
    for e in a[-1]:
        idx[e].append(i)

ans = 0
bs = list(map(lambda x: int(x) - 1, input().split()))
for b in bs:
    for id in idx[b]:
        cnt[id] -= 1
        if cnt[id] == 0:
            ans += 1
    print(ans)

