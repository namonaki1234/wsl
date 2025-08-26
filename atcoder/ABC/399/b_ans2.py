import io
import sys

# 下記に標準入力を記載
_INPUT = """\
4
3 12 9 9
"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

n = int(input())
p = list(map(int, input().split()))
for i in range(n):
  rank = 1
  for j in range(n):
    if p[j] > p[i]:
      rank += 1
  print(rank)
