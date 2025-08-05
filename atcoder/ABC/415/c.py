import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import combinations,permutations
import math

# 下記に標準入力を記載
_InPUT = """\
5
3
0010000
3
0010110
1
1
2
100
4
001110010101110

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

t = int(input())
for _ in range(t):
    # cnt = 0
    n = int(input())
    s = input().strip()

    # 1の位置を取得
    indices = [i+1 for i, char in enumerate(s) if char == '1']

    print(" ".join(map(str, indices)))

    # judege = True

    # for i in range(2**n):
    # for i in range(1<< n):

    num = (1 << n) - 1
    # for i in range(0,num+1,+1):
    #     print(format(i, '03b'))

    # for i in permutations(range(n), 2):
        # print(i)

    for order in permutations(range(n)):
        bits = [0] * n
        path = []

        for i in order:
            bits[i] = 1
            path.append(''.join(map(str, bits)))


        print(" -> ".join(path))
        if s not in path:
            print("Yes")
            break
    else:
        print("No")
