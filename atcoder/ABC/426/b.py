import io
import sys
from collections import defaultdict,deque,Counter

# 下記に標準入力を記載
_InPUT = """\
wwwwwwwwwv

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

s = input().strip()

counter_s = Counter(s)

# print(counter_s)
for k,v in counter_s.items():
    if v == 1:
        print(k)
