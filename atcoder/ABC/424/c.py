import io
import sys

# 下記に標準入力を記載
_InPUT = """\
6
6 1
2 3
3 1
0 0
4 6
0 0
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())

skill_set = set()
a_b = []
for i in range(n):
    a_b.append(list(map(int,input().split())))
    if a_b[-1][0] == 0 and a_b[-1][1] == 0:
        skill_set.add(i+1)
# print(skill_set)
for i in range(n):
    a, b = a_b[i]
    if a in skill_set or b in skill_set:
        skill_set.add(i+1)

print(len(skill_set))