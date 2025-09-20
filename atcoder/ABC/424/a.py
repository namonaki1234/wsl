import io
import sys
from collections import defaultdict,deque,Counter

# 下記に標準入力を記載
_InPUT = """\
10 10 10

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

a,b,c = map(int,input().split())

if a == b or b == c or c == a:
    print("Yes")
else:
    print("No")