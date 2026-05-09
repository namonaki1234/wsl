import io
import sys
from collections import Counter

# 下記に標準入力を記載
_InPUT = """\
1
5 100 200 300 400 500
1 5
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())

a = []
for i in range(n):
    row = [z for z in input().split()][1:]
    a.append(row)

x,y = map(int,input().split())

print(a[x-1][y-1])