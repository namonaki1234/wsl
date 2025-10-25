import io
import sys

# 下記に標準入力を記載
_InPUT = """\
8
1 (
2
1 (
1 )
2
1 (
1 )
1 )
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# sは最初空文字列
s = ""

q = int(input())


for _ in range(q):
    query = input().strip()
    if query[0] == '1':
        _, c = query.split()
        if c == ')':
            c = "1"
        else:
            c = "0"
        s += c
        # print(s)
    elif query[0] == '2':
        s = s[:-1]
    if s == "":
        print("Yes")
    else:
        index_0 = [i for i, ch in enumerate(s) if ch == '0']
        index_1 = [i for i, ch in enumerate(s) if ch == '1']
        for i in range(len(index_0)):
            if len(index_1) == 0 or index_1[0] == 0:
                print("No")
                break
            if len(index_1) == len(index_0) and index_0[i] < index_1[i]:
                print("Yes")
                break
        else:
            print("No")
