import io
import sys
from collections import Counter

# 下記に標準入力を記載
_InPUT = """\
8
72 74 69 70 73 75 71 77
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())
t = list(map(int, input().split()))

t_sorted = sorted(t)
# print(t_sorted)
cnt = 0
ans = []
for t_i in t_sorted:
    cnt += 1
    ans.append(t.index(t_i)+1)

    # if cnt == 3:
    #     exit()
ans = ans[:3]
print(*ans)
# t_index = []
# for i in range(n):
#     t_index.append(t[i])

# print(t_index)

# for t_item 
