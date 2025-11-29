import io
import sys
from collections import Counter

# 下記に標準入力を記載
_InPUT = """\
10 5
4 92
1 16
3 77
4 99
2 89
3 8
1 40
5 56
1 40
4 77
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

N, M = map(int, input().split())

# M 種類あるので、0-index にしたくない場合は M+1 個つくる
sumB = [0] * (M + 1)   # 各種類の大きさの合計
cnt  = [0] * (M + 1)   # 各種類の出現数

for _ in range(N):
    a, b = map(int, input().split())
    sumB[a] += b
    cnt[a]  += 1

# 平均値を出力
for k in range(1, M + 1):
    print(sumB[k] / cnt[k])
