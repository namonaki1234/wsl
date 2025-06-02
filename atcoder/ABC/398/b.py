import io
import sys
from collections import defaultdict,deque,Counter

# 下記に標準入力を記載
_INPUT = """\
2 2 4 4 3 3 3

"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

A = list(map(int, input().split()))

Counter_A = Counter(A)
x =[]
y =[]
for item in Counter_A.items():
    # print(item)
    if item[1] >= 3:
        x.append(item[0])
        
        
    elif item[1] == 2:
        y.append(item[0])
    
       

if len(x) >= 2:
            print("Yes")
            

elif len(x) == 1 and len(y) >= 1:
            print("Yes")

    
else:  
    print("No")  
    exit()  
# print(Counter_A)
# print(x)
# print(y)
# print(len(y))












    











