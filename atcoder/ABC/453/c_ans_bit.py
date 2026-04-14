import io
import sys
import itertools
# 下記に標準入力を記載
_InPUT = """\
5
2 5 2 2 1
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# 初期座標0.5
# i回目に正負どちらかにl_i移動
# 最大0を何回通れるか
# n<=20の全探索、分岐の数は2**n

n = int(input())
l = [int(l_i) for l_i in input().split( )]
# l = list(map(int, input().split()))

# 2倍して整数で扱う
l = [x * 2 for x in l]

ans = 0

# 2^n 通りの方向選択を全探索(0-2^(n-1))
for mask in range(1 << n):
    pos = 1   # 初期位置 0.5 を2倍したので 1
    cnt = 0   # この選び方で 0 を通り過ぎた回数

    for i in range(n):
        # maskは10進法でi回シフト演算をした後、→(mask >> i)で右からi桁目(0桁始まり)のbitを抜き出している(右に持ってきている)、その後 & 1とすることで抜き出したbitの一番右だけが取得できる、これはbigendian？
        if (mask >> i) & 1:
            # 正方向
            npos = pos + l[i]
        else:
            # 負方向
            npos = pos - l[i]

        # 0をまたいだか判定
        if pos * npos < 0:
            cnt += 1

        pos = npos

    ans = max(ans, cnt)

print(ans)