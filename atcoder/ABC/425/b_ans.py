import io
import sys
from itertools import permutations

# 下記に標準入力を記載
_InPUT = """\
7
3 -1 4 -1 5 -1 2

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載


n = int(input())
a = list(map(int, input().split()))
# for i in permutations([i + 1 for i in range(n)])で全探索を実装できる(1..7から7..1までの順列を全て試す)
for p in permutations([i + 1 for i in range(n)]):
    ok = True
    for i in range(n):
        #booleanには&=で代入できる
        # a[i] == -1 なら数字は何でもいいし、そうでなければ a[i] と p[i] が等しい必要があり、どちらかの条件を満たしていればいいのでorでつなぐ
        ok &= a[i] == -1 or p[i] == a[i]
    if ok:
        print("Yes")
        print(*p)
        break
else:
    print("No")
