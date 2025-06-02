import io
import sys
import numpy as np
from collections import deque

# 下記に標準入力を記載
_INPUT = """\
4
4 3 2 1
2 3 1 4


"""
sys.stdin = io.StringIO(_INPUT)


N=int(input())
P=[0]+list(map(int,input().split()))
Q=[0]+list(map(int,input().split()))

ans=[0]*(N+1)
for i in range(1,N+1):
  ans[Q[i]]=Q[P[i]]

print(*ans[1:])