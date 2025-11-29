import io
import sys


# 下記に標準入力を記載
_InPUT = """\
100 100
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# 高橋君の体重はw[kg]、風船をn個つけた物体の質量がnB[g]未満の時、空に飛び立つ
# その時最低何個

w,b = map(int,input().split())

ans = int((w *1000)/b +1)

print(ans)