import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np

# 下記に標準入力を記載
_InPUT = """\
10 22
47 81 82 95 117 146 165 209 212 215

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n,s = map(int, input().split())
t = list(map(int, input().split()))

sleep_time = s + 0.5
now_time = 0

for i in range(0,n):
    if i == 0:
        if t[i] - now_time < sleep_time:
            now_time = t[i]
        else:
            print("No")
            exit()
    if i == n-1:
        if t[n-1]-t[n-2] < sleep_time:
            now_time = t[i]
        else:
            print("No")
            exit()
    if i != n-1 and i != 0:
        if t[i+1]-t[i] < sleep_time:
            now_time = t[i+1]
        else:
            print("No")
            exit()

print("Yes")


    


    











