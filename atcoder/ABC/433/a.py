import io
import sys


# 下記に標準入力を記載
_InPUT = """\
1 100 2
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
#高橋ｘ歳、青木ｙ歳
#高橋の年齢がアオキの年齢のz倍になるかどうか判定
x,y,z = map(int, input().split())

a = (x-y*z)/(z-1)

if a >= 0 and a.is_integer():
    print("Yes")
else:
    print("No")
