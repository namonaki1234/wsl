import io
import sys
from collections import Counter

# 下記に標準入力を記載
_InPUT = """\
5 5
1 3 4 2 5
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n,m = map(int, input().split())
f = [int(x) for x in input().split(" ")]

#質問1
f_set =set(f)
if len(f) == len(f_set):
    print("Yes")
else:
    print("No")

#質問2
# f_counter = Counter(f)
# print(f_counter.values())
if len(f_set) == m:
    print("Yes")
else:
    print("No")