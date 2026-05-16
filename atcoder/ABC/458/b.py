import io
import sys
from collections import Counter

# 下記に標準入力を記載
_InPUT = """\
1 3
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

h,w = map(int,input().split())

if h == w == 1:
    print("0")
    exit()

if h == 1:
    start = "1" + "2" * (w-2) + "1"
    list_start =list(start)
    print(*list_start)
    exit()

if w == 1:
    start = "1"
    list_start =list(start)
    print(*list_start)
    for i in range(h-2):
        middle = "2"
        list_middle = list(middle)
        print(*middle)
    list_end = list_start
    print(*list_end)
    exit()

start = "2" + "3" * (w-2) + "2"
list_start =list(start)
print(*list_start)
for i in range(h-2):
    middle = "3"+"4"*(w-2)+"3"
    list_middle = list(middle)
    print(*middle)
list_end = list_start
print(*list_end)