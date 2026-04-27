import io
import sys
from collections import Counter

# 下記に標準入力を記載
_InPUT = """\
3 2
.#
#.
##
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

h,w = map(int, input().split())
# s = [input() for _ in range(h)]

# S = [
#     "",       # 0番目のダミー
#     " .#",   # 1行目
#     " #.",   # 2行目
#     " ##"    # 3行目
# ]
s = [""] # 縦の番兵
for i in range(h):
    s.append(" "+input()) # 横の番兵

ans = 0

for h1 in range(1,h+1):
    for h2 in range(h1,h+1):
        for w1 in range(1,w+1):
            for w2 in range(w1,w+1):
                ok = True
                #長方形内も点対称かチェック
                for i in range(h1,h2+1):
                    for j in range(w1,w2+1):
                        n_i = h1 + h2 - i
                        n_j = w1 + w2 - j
                        if s[i][j] != s[n_i][n_j]:
                            ok = False
                            break
                        
                if ok:
                    ans += 1

print(ans)