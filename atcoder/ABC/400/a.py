from collections import defaultdict
import io
import sys

# 下記に標準入力を記載
_INPUT = """\
400


"""


sys.stdin = io.StringIO(_INPUT)

A= int(input())

# div_400 =[1,2,4,5,8,10,16,20,25,40,50,80,100,200,400]

C = 400 % A
if C ==0 :
    print(int(400/A))
else:
    print(-1)


