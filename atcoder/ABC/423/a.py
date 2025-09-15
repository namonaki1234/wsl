import io
import sys
from collections import defaultdict,deque,Counter

# 下記に標準入力を記載
_InPUT = """\
10000000 24

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

x, c = map(int, input().split())

# for i in range(x, 0, -1000):
#     if i + (i // 1000) * c <= x and (i >= 1000):
#         print(i)
#         exit()
max_money = 0
for i in range(1000,x,1000):
    if i + (i // 1000) * c <= x and (i >= 1000):
        max_money = i
    # else:
    #     print(0)
    #     exit()
print(max_money)
