import io
import sys

# 下記に標準入力を記載
_INPUT = """\
4 1
1 2
2 3
2 4
1 2

"""

sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載



N,Q = map(int, input().split())
li = [[0] for _ in range(N)]

for i in range(N-1):
    a,b = map(int, input().split())
    li[a-1].append(b-1)
    li[b-1].append(a-1)
# print(li)

from collections import deque
d = deque([0])
dist = [-1] * N
dist[0] = 0

while d:
    x = d.popleft()
    for i in li[x]:
        # print(x,i,d)
        if dist[i] == -1:
            dist[i] = dist[x] + 1
            d.append(i)
# print(dist,li)
for _ in range(Q):
    a,b = map(int, input().split())
    a = dist[a-1]%2
    b = dist[b-1]%2
    if a!= b:
        print("Road")
    else:
        print("Town")



        
        
        


