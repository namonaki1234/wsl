import io
import sys
from collections import defaultdict,deque,Counter

# 下記に標準入力を記載
_INPUT = """\
7 3

"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載


N,M= map(int, input().split())

X=(N**(M+1)-1)/(N-1)
if X <= 10**9:
    print(int(X))
else:
    print("inf")


    











