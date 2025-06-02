import io
import sys

# 下記に標準入力を記載
_INPUT = """\
7
13 17 14 20 18 27 10
"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

N = int(input())
d = map(int,input().split())
result = sum(d)/N
result_str = str(result) 
result_int = int(result)
if int(result_str[result_str.index('.') + 1:][-1]) == 0:    
        print(result_int)
else:
    print(result_int+1)