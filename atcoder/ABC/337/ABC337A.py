
n = int(input())

t = 0 # チーム高橋の合計得点
a = 0 # チーム青木の合計得点

for i in range(n):
    x,y= map(int, input().split())
    
    # それぞれのチームに得点をプラスする
    t += x
    a += y

# 合計得点の大小に従って結果を出力
if t > a:
    print('Takahashi')
elif t < a:
    print('Aoki')
else:
    print('Draw')




'''
import numpy as np
import re 
n = int(input())

a = []
i = 0
while i < n :
    a.append(list(map(int,input().split())))

    i += 1

for u in a:
    print(u)
    
print(a)

b = []  
c = []

k = 0
while k < n-1:
    b.append(a[k][0]+ a[k+1][0])
    c.append(a[k][1]+ a[k+1][1])
    k += 1
    
print(b)
print(c)

'''