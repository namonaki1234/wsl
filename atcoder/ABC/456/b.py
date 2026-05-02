import io
import sys
from collections import Counter

# 下記に標準入力を記載
_InPUT = """\
1 2 3 4 5 6
1 2 3 4 5 6
1 2 3 4 5 6
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

a = [[0] * 6 for _ in range(3)]
visited = [[False] * 6 for _ in range(3)]

for i in range(3):
    a[i] = input().split(" ")


con1 = a[0].count("4") and a[1].count("5") and a[2].count("6")

con2 = a[0].count("4") and a[1].count("6") and a[2].count("5")
con3 = a[0].count("5") and a[1].count("4") and a[2].count("6")
con4 = a[0].count("5") and a[1].count("6") and a[2].count("4")
con5 = a[0].count("6") and a[1].count("4") and a[2].count("5")
con6 = a[0].count("6") and a[1].count("5") and a[2].count("4")

ans = 0.0
if con1:
    ans += (a[0].count("4") * a[1].count("5") * a[2].count("6")) / 216
    # print((a[0].count("4") * a[1].count("5") * a[2].count("6")) / 216)
if con2:
    ans += (a[0].count("4") * a[1].count("6") * a[2].count("5")) / 216
    # print((a[0].count("4") * a[1].count("6") * a[2].count("5")) / 216)
if con3:
    ans += (a[0].count("5") * a[1].count("4") * a[2].count("6")) / 216
if con4:
    ans += (a[0].count("5") * a[1].count("6") * a[2].count("4")) / 216
if con5:
    ans += (a[0].count("6") * a[1].count("4") * a[2].count("5")) / 216
if con6:
    ans += (a[0].count("6") * a[1].count("5") * a[2].count("4")) / 216
    
print(ans)
# print(con2)
# for i in range(1,4):
#     for j in range(1,7):
# print(a)
# print(a[0].count("4"))
