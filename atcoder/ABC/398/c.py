import io
import sys
from collections import defaultdict,deque,Counter

# from more_itertools import chunked

# 下記に標準入力を記載
_INPUT = """\
9
2 9 9 7 9 2 4 5 8

"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載
N = int(input())
A =list(map(int, input().split()))
Counter_A = Counter(A)
# print(Counter_A)
num_max = 0
unique_num = []
for item in Counter_A.items():
    # print(item)
    # ユニーク
    if item[1] == 1:
        unique_num.append(item[0])
        num_max = max(num_max, unique_num[-1])
    
# print(num_max)
for i in range(N):
    if A[i] == num_max:
        print(i+1)
        break
    else:
        continue
# print(unique_num)
if len(unique_num) == 0:
    print(-1)
    exit()











