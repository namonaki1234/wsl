import io
import sys

# 下記に標準入力を記載
_INPUT = """\
5
3 9 5 3 1

"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

n = int(input())
a = list(map(int, input().split()))

ans = n + 1
# 1_000_000の_は桁区切りで百万には変わりない
#　aの値をpositionつまり位置で管理する
pos = [[] for _ in range(1_000_001)]

for i in range(n):
    pos[a[i]].append(i)

for i in range(1_000_001):
    # どこかの数字がダブるまではmin関数に入れる必要がないのでスキップする
    if len(pos[i]) < 2:
        continue
    for j in range(len(pos[i]) - 1):
        ans = min(ans, pos[i][j + 1] - pos[i][j] + 1)

if ans == n + 1:
    print(-1)
else:
    print(ans)
