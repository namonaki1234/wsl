import io
import sys

# 下記に標準入力を記載
_InPUT = """\
5 5
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

m,d = map(int,input().split())

goseku = [(1,7),(3,3),(5,5),(7,7),(9,9)]   

if (m,d) in goseku:
    print("Yes")
else:
    print("No")

