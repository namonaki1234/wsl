import io
import sys
import numpy as np
import math
import itertools

# 下記に標準入力を記載
_INPUT = """\
3
3 1 2 3
4 1 2 2 1
6 1 2 3 4 5 6

"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

N = int(input())
K = [0] * N

A_list =[]
A =[]
for i in range(N):
    #上から順番に代入していく
    K[i], A = input().split(' ', 1)
    K[i] = int(K[i])
    
    A_list.append(list(map(int, A.split())))


n = 5
r = 3
combination = math.comb(n, r)

    
elements = range(1,len(K)+1)
r = 2

# 組み合わせを生成
combinations = list(itertools.combinations(elements, r)) 

for i in range(len(combinations)):
    combinations[i]









