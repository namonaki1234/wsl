import io
import sys
from collections import defaultdict

# 下記に標準入力を記載
_INPUT = """\
5
3 9 5 3 1

"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載
# 尺取り法（two pointers）

N = int(input())
A = list(map(int, input().split()))

cnt = defaultdict(int)
dup = 0          # 窓内で 2回以上出ている値の種類数
ans = N + 1

l = 0
for r, x in enumerate(A):
    cnt[x] += 1
    if cnt[x] == 2:
        dup += 1

    # 窓内に重複がある限り、左を詰めて最短化を狙う
    while dup > 0:
        ans = min(ans, r - l + 1)
        y = A[l]
        cnt[y] -= 1
        if cnt[y] == 1:   # 2→1 に落ちたら、重複の種類が1つ消える
            dup -= 1
        l += 1

print(-1 if ans == N + 1 else ans)
