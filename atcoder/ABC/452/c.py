import io
import sys
import math
from collections import Counter

# 下記に標準入力を記載
_InPUT = """\
5
5 3
5 2
4 1
5 1
3 2
8
retro
chris
itchy
tuna
crab
rock
cod
ash
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# n本の肋骨、1本の脊椎からなる魚の骨のオブジェ
# 肋骨iに長さA_iの文字列を刻む、B_i文字目は脊椎のi文字目と一致

n = int(input())
ab = [list(map(int, input().split())) for _ in range(n)]
m = int(input())
s = [input() for _ in range(m)]

def s_and_s_len_return(s,m):
    s_and_s_len = []
    for i in range(m):
        s_and_s_len.append((s[i],len(s[i])))
    return s_and_s_len
s_and_s_len = s_and_s_len_return(s,m)
print(s_and_s_len)
spine = ""
print(spine)

ans_criterion = [ans_i for ans_i in s if len(ans_i) == n]
# print(ans_criterion)
for i in range(len(ans_criterion)):
    a = s[i][ab[i][0]]
    b = s[i][ab[i][-1]]
    target = ans_criterion[i]
    a_criterion = [a_i for a_i in s if len(a_i) == a]
    count = 0
    for j in range(len(a_criterion)):
        if a_criterion[j][b] == target[i]:
            count += 1
        else:
            break

# s_counter = Counter(s)
# print(s_counter)

# print(n)
# print(ab)
# print(m)
# print(s)

# for i in range(n):
#     target = s[i]
#     for j in range(n):
#         b = s[j][ab[i][-1]]
#         if b == s[j][]
# for i in range(n):
#     spine += s[i][ab[i][-1]]

# s_copy = s.copy()
# visited = []
# end = 0
# k = 0
# for i in range(math.factorial(5)):
              

