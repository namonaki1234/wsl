import io
import sys

# 下記に標準入力を記載
_InPUT = """\
4
4 3 2 5
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())
a = list(map(int, input().split()))

a_copy = a.copy()
a.reverse()
# print(a)

iteration = n
cnt = 0
ans = []
start = 1
for i in range(n):
    for j in range(start,len(a)):
        if a[i] < a[j]:
            # print(a_copy.index(a[j]) + 1)
            ans.append(a_copy.index(a[j]) + 1)
            cnt += 1
            break
        if a[i] >= a[j]:
            ans.append(-1)
            break

        if cnt == len(a)-1-start+1:
            # print(-1)
            ans.append(-1)
    cnt = 0
    start += 1
    if start == len(a):
        ans.append(-1)
        break
for k in range(len(ans)-1,-1,-1):
    print(ans[k])
# max_save = a[0]
# second_save = -1
# for i in range(n):
#     max_save = max(max_save, a[i])
#     if a[i] >= max_save and a[i] >= second_save :
#         print(-1)
#     else:
#         if max_save > a[i] and second_save > a[i]:
#             print(a.index(second_save) + 1)
#         if max_save > a[i] and second_save < a[i]:
#             print(a.index(max_save) + 1)
#     second_save = max(second_save, a[i] if a[i] < max_save else -1)

