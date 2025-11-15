import io
import sys

# 下記に標準入力を記載
_InPUT = """\
903
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

x = int(input())

list_x = list(str(x))
list_x.sort()
# print(list_x)
if list_x[0] == '0':
    for i in range(1, len(list_x)):
        if list_x[i] != '0':
            list_x[0], list_x[i] = list_x[i], list_x[0]
            break
print(''.join(list_x))