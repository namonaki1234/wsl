import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np

# 下記に標準入力を記載
_InPUT = """\
3 1
100 160 130
120

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載


# 入力の受け取り
N,Q=map(int,input().split())
A=list(map(int,input().split()))

# Aをソート
A.sort()

# Q回
for i in range(Q):
    # 入力の受け取り
    x=int(input())

    # xがAの中で最小の要素以下である場合
    if x<=A[0]:
        # N人
        print(N)

    # xがAの中で最大の要素より大きい場合
    elif A[N-1]<x:
        # 0人
        print(0)

    # 二分探索
    else:
        # 左
        l=0
        # 右
        r=N-1

        # 1<右-左の間
        while 1<r-l:
            # 中央
            c=(l+r)//2

            # A[c]<xならば(条件を満たさない場合:x以上の排反つまり余事象)
            if A[c]<x:
                # 左=中央と更新
                l=c
            # そうでなければ(x≤A[c] 条件を満たす場合)
            else:
                # 右=中央と更新
                r=c

        # 答え：(N-r)人
        print(N-r)
