import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np

# 下記に標準入力を記載
_InPUT = """\
97 210


"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

l, r = map(int, input().split())

n = r - l + 1
cnt = 0
for num in range(l,r+1):

    # 二分探索
    num_str = str(num)
    max_num = 0
    for num_i in num_str:
        max_num = int(num_i,max_num)
        # 左
        l=0
        # 右
        r=n-1

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

print(cnt)