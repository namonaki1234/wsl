import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np

# 下記に標準入力を記載
_InPUT = """\
6
1 6
2 5
3 4
4 3
5 2
6 1


"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())

cnt = 0

for _ in range(n):
    a,b = map(int, input().split())
    if a < b:
        cnt += 1
print(cnt)