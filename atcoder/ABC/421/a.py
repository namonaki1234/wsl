import io
import sys
from collections import defaultdict,deque,Counter

# 下記に標準入力を記載
_InPUT = """\
2
wang
li
2 wang

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())
s = [input() for _ in range(n)]
x, y = input().split()
x = int(x)

if s[x-1] == y:
    print("Yes")
else:
    print("No")