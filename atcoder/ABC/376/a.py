import io
import sys

# 下記に標準入力を記載
_InPUT = """\
6 5
1 3 7 8 10 12
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n,c = map(int, input().split())
t = list(map(int, input().split()))

candy_count = 0
candy_get_time = []
for i in range(n-1):
    if t[i+1] - t[i] < c:
        if i == 0:
            candy_count += 1
        continue
    elif i == 0 or t[i] - t[i-1] >= c:
        candy_count += 1

print(candy_count)