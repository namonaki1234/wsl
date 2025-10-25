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

def check_0_1(s):
    if len(s) % 2 != 0 or s == "":
        return "No"
    count_0 = 0
    count_1 = 0
    count_chunk = 0
    for c in s:
        # 0のとき
        if c == '0':
            count_0 += 1
        # 1のとき
        else:
            if count_0 == 0:
                count_1 += 1
            count_0 -= 1
            count_chunk += 1
    if len(s)/2 == count_chunk and count_0 == 0 and count_1 == 0:
        return "Yes"
    else:
        return "No"

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
        print(check_0_1(s))
