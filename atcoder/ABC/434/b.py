import io
import sys
from collections import Counter

# 下記に標準入力を記載
_InPUT = """\
10 5
4 92
1 16
3 77
4 99
2 89
3 8
1 40
5 56
1 40
4 77
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n,m = map(int, input().split())

a=[]
b=[]
ab_dict = {}
ab_list =[]
for i in range(n):
    a_i,b_i = map(int, input().split())
    # ab_dict.append(a_i,b_i)
    # ab_dict[a_i]=b_i
    # ab_list.append({a_i:b_i})
    ab_list.append((a_i,b_i))
    # a.append(a_i)
    # b.append(b_i)

# print(ab_dict)
# print(ab_list)
# print(a,b)

ab_list.sort()

# print(ab_list)

cnt = 1
sum = 0
avg_list = []
avg_cnt = 0

for ab in ab_list:
    if cnt != ab[0]:
        cnt+=1
        avg_list.append(sum/avg_cnt)
        sum = 0
        avg_cnt= 0
    if cnt == ab[0]:
        sum += ab[1]
        avg_cnt+=1
        if cnt == m:
            avg_list.append(sum/avg_cnt)

for i in range(m):
    print(avg_list[i])

# i=0
# while i<n:
#     for ab in ab_list:
#         if cnt == ab[0]:
#             sum += ab[1]
#             avg_cnt+=1
#             continue
#     avg_list.append(sum/avg_cnt)
#     cnt +=1
#     sum = 0
#     i+=1

# ab_list_cnt = Counter(ab_list)
# print(ab_list_cnt)

