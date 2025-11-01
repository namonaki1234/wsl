import io
import sys
from collections import Counter
from itertools import combinations
import math

# 下記に標準入力を記載
_InPUT = """\
11 4 2
abbaaabaaba
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n, a, b = map(int, input().split())
s = input().strip()

a_index = []
b_index = []
for i in range(n):
    if s[i] == 'a':
        a_index.append(i)
    elif s[i] == 'b':
        b_index.append(i)
print(a_index)
print(b_index)

# s_slice = s[4:9]
# print(s_slice)
a_slice = []
b_slice = []
for i in range(0,n-(a-1)):
    if len(a_index[i:i+a]) < a:
        break
    a_slice.append(a_index[i:i+a])

b_copy = b
b -= 1
while b > 0:
    for j in range(0,n-(b-1)):
        if len(b_index[j:j+b]) < b:
            break
        b_slice.append(b_index[j:j+b])
    b -= 1
print(a_slice)
print(b_slice)

count = 0
ans_count = 0
for a_s in a_slice:
    count = 0
    for b_s in b_slice:
        if a_s[0] < b_s[0] and a_s[-1] > b_s[-1]:
            # if  count < b_copy - 1:
            count += 1
            print(b_s[0], b_s[-1])
            # else:
            #     continue
    print(count)
    if count < b_copy:
        ans_count += 1
print(ans_count)