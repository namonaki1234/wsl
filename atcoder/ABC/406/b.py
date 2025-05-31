import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np

# 下記に標準入力を記載
_InPUT = """\
8
19 5 5 19 5 19 4 19


"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

m = int(input())
c = list(map(int, input().split())) 

set_c = set(c)

sorted_c = sorted(set_c)

print(len(sorted_c))
# print(" ".join(map(str, sorted_c)))
print(*sorted_c)   


    











