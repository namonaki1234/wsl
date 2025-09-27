import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import combinations

# 下記に標準入力を記載
_InPUT = """\
8 7
.######
##....#
#.###.#
#.#.#.#
#.#.#.#
#.#####
#...#..
#####..
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
#.が白、#が黒
# すべての黒マスについて、上下左右の黒マスの個数が 2,4 のいずれかであるかを調べる

H, W = map(int, input().split())

# 上下左右に 1 マスずつ白マス(.)を広げたマス目を用意する(番兵（sentinel）)←これにより、境界のマスの処理がエラーにならないようにする
# 「番兵（sentinel）」というのは、本来のデータの外側にあえて特別な値を置いておくことで、境界条件の分岐処理を単純化するテクニックのことです。
# 由来（言葉の意味）
# sentinel は英語で「見張り」「番兵」という意味。
# 配列やグリッドの外に置いておく「見張り役のデータ」が、外に出ようとしたときに必ずキャッチしてくれる。だから番兵と呼ぶ。
field = ['.' * (W + 2)] + ['.' + input() + '.' for i in range(H)] + ['.' * (W + 2)]

for i in range(1, H + 1):
    for j in range(1, W + 1):
        black_count = 0 # 周囲の黒マスの個数
        if field[i - 1][j] == '#':
            black_count += 1
        if field[i][j - 1] == '#':
            black_count += 1
        if field[i + 1][j] == '#':
            black_count += 1
        if field[i][j + 1] == '#':
            black_count += 1

        if field[i][j] == '#' and black_count != 2 and black_count != 4:
            # 条件を満たさないマスがあったら、No を出力して終了
            print('No')
            exit(0)

# すべてのマスが条件を満たすので、Yes を出力して終了
print('Yes')
