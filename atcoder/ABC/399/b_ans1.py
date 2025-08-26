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
rank = [0 for i in range(n)]
r = 1
while True:
  x = 0
  for i in range(n):
    if rank[i] != 0:
      continue
    x = max(x, p[i])
  if x == 0:
    break
  k = 0
  for i in range(n):
    if p[i] == x:
      rank[i] = r
      k += 1
  r += k

for i in range(n):
  print(rank[i])
