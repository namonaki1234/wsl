import io
import sys
from collections import defaultdict,deque,Counter

# 下記に標準入力を記載
_InPUT = """\
10

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())
ans = 0
for i in range(1,n+1):
    ans += (-1)**i * (i ** 3)
print(ans)