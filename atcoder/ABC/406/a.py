import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np

# 下記に標準入力を記載
_InPUT = """\
22 40 22 30

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

a,b,c,d = map(int, input().split())

deadline = a*60+b
submit = c*60+d
if deadline >= submit:
    print("Yes")
else:
    print("No")


    


    











