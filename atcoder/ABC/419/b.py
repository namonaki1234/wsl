import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import combinations

# 下記に標準入力を記載
_InPUT = """\
8
1 5
1 1
1 1
1 9
2
2
1 2
2


"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

q = int(input())

ball_bag_cnt =[]
for _ in range(q):
    string = list(map(int, input().split()))
    length = len(string)
    # タイプ1
    if length == 2:
        type,x = string
        ball_bag_cnt.append(x)
        # print(ball_bag_cnt)
    # タイプ2
    elif length == 1:
        type = string[0]
        min_ball = min(ball_bag_cnt) if ball_bag_cnt else 0
        print(min_ball)
        ball_bag_cnt.remove(min_ball)
