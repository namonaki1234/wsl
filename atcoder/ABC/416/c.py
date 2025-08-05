import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import combinations,permutations
import math

# 下記に標準入力を記載
_InPUT = """\
3 2 6
abc
xxx
abc


"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n, k, x = map(int, input().split())

s = []
for _ in range(n):
    s.append(input().strip())


string_s = []
for i in range(1,n+1):
    for j in range(1,n+1):
        string_s.append(s[i-1] + s[j-1])

string_s = sorted(string_s)

print(string_s[x-1])
