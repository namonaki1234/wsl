import io
import sys


# 下記に標準入力を記載
_InPUT = """\
10 5
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# 一年d日、その年の最初のコンテストがf日目

d,f = map(int,input().split())

cnt_con_res = (d-f) % 7
ans =  7 - cnt_con_res

print(ans)
