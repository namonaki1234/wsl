import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import islice

# 下記に標準入力を記載
_InPUT = """\
12
11 3 8 9 5 2 10 4 1 6 12 7

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())
p = list(map(int, input().split())) 

# combinations = list(itertools.combinations(p, 3))
# print(combinations)
cnt_max = False
cnt_min = False

arr_3 = [list(islice(p, i, i + 3)) for i in range(len(p) - 2)]
# print(arr_3)

cnt_childer = 0
for i in range(len(arr_3)):
        if arr_3[i][0] < arr_3[i][1] > arr_3[i][2]:
            cnt_max = True
            max_head = arr_3[i][0]
        if arr_3[i][0] > arr_3[i][1] < arr_3[i][2]:
            cnt_min = True
            min_head = arr_3[i][0]
        if cnt_max == True and cnt_min == True:
            if  p.index(max_head) < p.index(min_head):
                if p[p.index(max_head)-1] < max_head:
                    cnt_childer += 2
                else :cnt_childer += 1
            else:
                if p[p.index(min_head)-1] < min_head and p.index(min_head)-1 > 0:
                    cnt_childer += 1
                cnt_childer += 0

print(cnt_childer)

    











