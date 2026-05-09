import io
import sys

# 下記に標準入力を記載
_INPUT = """\
8
5 7 5 7 7 5 7 7
"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

n = int(input())
h = [x for x in input().split(" ")]

dp = [[1]*n for _ in range(n)]

ans = 1

for i in range(n):
     # i + d が N 未満になる範囲だけ見る
    for d in range(1,n-i):
        next_i = i + d

        if h[next_i] == h[i]:
            dp[next_i][d] = dp[i][d] + 1
            ans = max(ans,dp[next_i][d])

print(ans)