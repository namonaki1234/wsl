import io
import sys
from collections import Counter

# 下記に標準入力を記載
_InPUT = """\
6 5
1
3
2
3 5
5
5 3 1 4 2
5
5 1 3 4 2
5
3 4 1 5 2
5
5 1 3 2 4
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# lはリストの長さ（選べるジュースの数）、xはジュースのid識別番号

n,m = map(int,input().split())

juice_list = [i for i in range(1,m+1)]
# print(juice_list)
for i in range(n):
    l = int(input())
    x = list(map(int, input().split()))
    juice_set = set(juice_list)

    x_and_list = [v for v in x if v in juice_set]
    
    # print(juice_list.pop(0))
    # print(x_and_list.pop(0))
    if len(x_and_list) == 0:
        print(0)
        continue
    print(x_and_list[0])
    juice_list.remove(x_and_list[0])

