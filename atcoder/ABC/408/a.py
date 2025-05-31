import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np

# 下記に標準入力を記載
_InPUT = """\
22 11


"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

a,b = map(int, input().split())

standard = int((a/b))
diff = abs(standard-a/b)

if diff >= 0.5:
    print(int(standard + 1))
elif diff < 0.5:
    print(int(standard))
elif a % b == 0:
    if a < b:
        print(b/a)
    else:
        print(a/b)
