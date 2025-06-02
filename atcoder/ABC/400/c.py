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
    for b in range(1, int(math.sqrt(X))+1):
        if  X == 2**(math.log2(X/b**2)) *b**2 and math.log2(X/b**2).is_integer():
            ans.append((X))
            break
        
print(len(ans))








