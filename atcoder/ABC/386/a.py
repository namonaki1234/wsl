import io
import sys
from collections import defaultdict,deque,Counter

# 下記に標準入力を記載
_InPUT = """\
1 3 4 1

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

a, b, c, d = map(int, input().split())

# print(Counter([a, b, c, d]))

counter = Counter([a, b, c, d])

if len(counter) != 2:
    print("No")
    exit()

# print(counter.values())
# print(counter.keys())

if (list(counter.values())[0] == 3 and list(counter.values())[1] == 1) or (list(counter.values())[1] == 3 and list(counter.values())[0] == 1) \
    or (list(counter.values())[0] == 2 and list(counter.values())[1] == 2):
    print("Yes")
else:
    print("No")