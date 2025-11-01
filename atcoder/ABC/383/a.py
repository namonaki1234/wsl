import io
import sys


# 下記に標準入力を記載
_InPUT = """\
3
1 8
10 11
21 5
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())

volume = 0
t = []
v = []
for _ in range(n):
    t_i, v_i = map(int, input().split())
    t.append(t_i)
    v.append(v_i)
# print(t, v)

for i, t_i in enumerate(t):
    if t_i == 1:
        # print(volume)
        volume += v[0]
    else:
        if i == 0:
            volume += v[0]
            continue
        else:
            volume += - (t[i] - t[i-1]) * 1
        if volume < 0:
                volume = 0
        volume += v[i]
        # print(volume)

print(volume)