import io
import sys
from collections import defaultdict,deque,Counter

# 下記に標準入力を記載
_InPUT = """\
5 8
1 2
1 3
1 4
2 3
2 5
3 4
3 5
4 5


"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
sys.setrecursionlimit(10**9)

N, M = map(int, input().split())

edges = [tuple(map(int, input().split())) for _ in range(M)]

ans = M
# 2^N 通りの塗り方を全部探索する
for bit in range(2 ** N):
    delete_count = 0
    for u, v in edges: # それぞれの辺を見て
        if (1 & (bit >> u)) == (1 & (bit >> v)): # 結んでいる頂点が同じ色で塗られていたら
                delete_count += 1 # カウントを増やす
    ans = min(ans, delete_count)

print(ans)
