import io
import sys
from collections import defaultdict,deque,Counter

# 下記に標準入力を記載
_INPUT = """\
6
1 30
1 1
1 100
2
1 3
2

"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

N = int(input())
Q = [input()for _ in range(N)]

human =deque()
num =deque()
for i in range(N):
    if Q[i][0] == "1":
        human.append(1)
        if len(Q[i]) == 5:
            num.append(int(Q[i][2:5]))
        elif len(Q[i]) == 4:
            num.append(int(Q[i][2:4]))
        elif len(Q[i]) == 3:
            num.append(int(Q[i][2:3]))
        else:
            num.append(int(Q[i][2]))
    else:
        human.popleft()
        print(num.popleft())
        




    


    











