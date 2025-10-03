import io
import sys

# 下記に標準入力を記載
_InPUT = """\
7
3 -1 4 -1 5 -1 2

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())
a = list(map(int,input().split()))

idx_not_minus = []
not_minus_num = []
for i in range(n):
    if a[i] != -1:
        idx_not_minus.append(i)
        if a[i] in not_minus_num:
            print("No")
            exit()
        not_minus_num.append(a[i])

# print((idx_not_minus))
# print((not_minus_num))
p = [-1] * n
range_n = [i for i in range(1,n+1) if i not in not_minus_num]
for i in range(1,n+1):
    if i in not_minus_num:
        p[i-1] = a[idx_not_minus[0]]
        idx_not_minus.pop(0)
for p_i in p:
    if p_i == -1:
        p_i = range_n.pop()
        p[p.index(-1)] = p_i
        # p[i] = a[idx_not_minus[0]]
        # ,p[a[idx_not_minus[0]]-1] = , p[i]

if len(set(p)) != n:
    print("No")
    exit()
print("Yes")
print(*p)