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

H, W = map(int, input().split())

S = [""]  # ダミー
for _ in range(H):
    S.append(" " + input())  # 列も1-indexedにするため，先頭に空白を足す

ans = 0

for h1 in range(1, H + 1):
    for h2 in range(h1, H + 1):
        for w1 in range(1, W + 1):
            for w2 in range(w1, W + 1):

                ok = True

                for i in range(h1, h2 + 1):
                    for j in range(w1, w2 + 1):
                        ni = h1 + h2 - i
                        nj = w1 + w2 - j

                        if S[i][j] != S[ni][nj]:
                            ok = False

                if ok:
                    ans += 1

print(ans)