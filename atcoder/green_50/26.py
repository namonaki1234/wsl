import io
import sys

# 下記に標準入力を記載
_INPUT = """\
a
4
2 1 p
1
2 2 c
1

"""

sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

#・文字列操作の計算量を理解する ・dequeの使い方を身につける
S=input()
Q=int(input())

inv=False
 
from collections import deque
S_deque=deque()
 
for i in range(len(S)):
    S_deque.append(S[i])
 
for i in range(Q):
    TFC=list(map(str, input().split()))
 
    if TFC[0]=="1":
        if inv==False:
            inv=True
        else:
            inv=False

    else:
        F=TFC[1]
        C=TFC[2]

        if inv==False:
            if F=="1":
                S_deque.appendleft(C)
            else:
                S_deque.append(C)
        else:
            if F=="1":
                S_deque.append(C)
            else:
                S_deque.appendleft(C)
 
ans="".join(S_deque)
 
if inv==True:
    ans=ans[::-1]
 
print(ans)
