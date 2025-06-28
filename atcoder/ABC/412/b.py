import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np

# 下記に標準入力を記載
_InPUT = """\
aBCdE
abcdcba


"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

s = input().strip()
t = input().strip()

upper_cnt = 0

for i in range(len(s)):
    if s[i].isupper():
        upper_cnt += 1
        if upper_cnt == 2:
            if s[i-1] in t:
                print("Yes")
                exit()
if upper_cnt == 0:
    print("Yes")
    exit()
if upper_cnt == 1:
    if not s[0].isupper() and s[0] in t:
        print("Yes")
        exit()
print("No")
