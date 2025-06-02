import io
import sys
from collections import defaultdict,deque,Counter
import math
# from more_itertools import chunked

# 下記に標準入力を記載
_INPUT = """\
20

"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載
N = int(input())

ans =[]
for X in range(2, N+1,2):
        for i in range(1, int(math.sqrt(N))+1):
            if X == 2 and N >= 2:
                ans.append((X))
                break
            if (X/2)%(i**2)==0 and (X/2)/(i**2)%2==0 and (X)/(i**2) >= 2:
                ans.append((X))
                break
ans_set = set(ans)
print(len(ans_set))








