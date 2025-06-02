import io
import sys
from collections import defaultdict,deque,Counter

# 下記に標準入力を記載
_INPUT = """\
4
3 12 9 9
"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載


N= int(input())
P = list(map(int, input().split()))
index_P = list(enumerate(P))
# print(index_P)
ans = [0]*N
sorted_P = sorted(index_P,key=lambda x: x[1],reverse=True)
# print(sorted_P)
r = 1
indices, values = zip(*sorted_P)
print(indices)
print(values)
    
Counter_P = Counter(P)
print(Counter_P)
k =[]
for item in Counter_P.items():
    # print(item)
    # 出現回数がitem[1]の時
    k.append(item[1])
print(k)

for item in Counter_P.items():
    print(item)
    # 出現回数がitem[1]の時
    # k = item[1]
    print(k)
    for i in range(N):
        if item[0] == values[i]:
            if i== 0:
                ans[indices[i]] = r
            elif values[i] == values[i-1]:
                ans[indices[i]] = r
            else:
                        r += k[indices[i-1]]
                        ans[indices[i]] = r
                

print(*ans)
    
# set_P = set(P)
# print(set_P)
# ans = [0]*N
# for i in range(N):
#     if i == 0:
#         ans[i] = Counter_P[P[i]]
#     else:
#         ans[i] = P.count(set_P[i])+ans[i-1]
# print(ans)

# for item in Counter_P.items():
    
#     print(item)
#     if item[1] == 1:
#         r[item[1]] = 1
        
   
   

# print(max(Counter_P))
# for i,p in enumerate(P):
#     if p == max(P):
#         r[i] = 1
#         # r += Counter_P[1]
#     # print(i)
#     # print(p)


# print(r)






    











