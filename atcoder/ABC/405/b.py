import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np

# 下記に標準入力を記載
_InPUT = """\
10 4
1 3 3 4 2 1 3 1 2 4

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n,m = map(int, input().split())
a = list(map(int, input().split())) 

# a_sorted = sorted(a,reverse=True)
set_cnt_m = set()
index = 0
for i in range(1,m+1):
    if i not in a:
        print(0)
        exit()

for i,item in enumerate(a):
    set_cnt_m.add(item)
    if len(set_cnt_m) == m:
        index = i
        break

print(len(a)-(index))

    


    











