import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np

# 下記に標準入力を記載
_InPUT = """\
6
2 3 5 7 11 13
1


"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())
a = list(map(int, input().split()))
x = int(input())

if x in a:
    print("Yes")
else:
    print("No")
