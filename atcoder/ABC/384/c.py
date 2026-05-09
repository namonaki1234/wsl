import io
import sys


# 下記に標準入力を記載
_InPUT = """\
400 500 600 700 800
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

# a,b,c,d,e = map(int,input().split())
# keys = [a,b,c,d,e]
keys = [x for x in range(97,102)]
values = list(map(int,input().split()))
# print(a,b,c,d,e)

# ans_dict = dict.fromkeys(keys,abcde_list)
ans_dict = dict(zip(keys,values))
print(ans_dict)


