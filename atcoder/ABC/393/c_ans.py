from collections import defaultdict
import io
import sys

# 下記に標準入力を記載
_INPUT = """\
3 5
1 2
2 3
3 2
3 1
1 1

"""

sys.stdin = io.StringIO(_INPUT)

# NとMを取得
N, M = list(map(int, input().split()))
ans = 0
edge_count = defaultdict(int)

# 辺の入力と処理
for _ in range(M):
    u, v = map(int, input().split())
    
    # 自己ループのカウント
    if u == v:
        ans += 1
        continue
    
    # 順序を入れ替え
    if u > v:
        u, v = v, u
    
    # 辺の出現回数をカウント
    edge_count[(u, v)] += 1

# 重複辺のカウント
for k in edge_count.values():
    ans += k - 1

# 結果の出力
print(ans)
