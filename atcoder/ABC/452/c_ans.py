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

# 意味は：
# lpc[長さ][位置][文字番号]
# 長さは 1〜10
# 位置は 0〜9
# 文字番号は a→0, b→1, ..., z→25

lpc = [[[False] * 26 for _ in range(10)] for _ in range(11)]
# lpc = [length][position][char]

for i in range(m):
    for j in range(len(s[i])):
        ch = ord(s[i][j]) - ord('a')
        lpc[len(s[i])][j][ch] = True

# print(lpc)

ans = []

for i in range(m):
    spine = s[i]

    if len(spine) != n:
        ans.append("No")
        continue

    ok = True

    for j in range(n):
        a = ab[j][0]
        b = ab[j][-1] - 1
        c = ord(spine[j]) - ord("a") #aの文字コードが97なので、a-aつまり97-97=0としてa:0,b:1...のように変換している

        # 出力が一つでも存在するかという指定だから、余事象的に書くことでシンプルにしている
        # 一つでもfalseならok=Falseとなる、つまりすべてのlpc[a][b][c]がTrueでないといけない
        if not lpc[a][b][c]:
            ok = False
            break
    
    if ok:
        ans.append("Yes")
    else:
        ans.append("No")

print(*ans, sep="\n")