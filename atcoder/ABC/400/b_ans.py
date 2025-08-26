import io
import sys
from collections import defaultdict,deque,Counter

# 下記に標準入力を記載
_INPUT = """\
998244353 99
"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載


N,M= map(int, input().split())
# inf は 「無限大 (infinity)

# X=(N**(M+1)-1)/(N-1)
# if X <= 10**9:
#     print(int(X))
# else:
#     print("inf")

x = sum(N**i for i in range(M+1))
print(x if x <= 10**9 else "inf")
