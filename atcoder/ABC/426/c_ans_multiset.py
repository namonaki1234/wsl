import io
import sys

# 下記に標準入力を記載
_InPUT = """\
8 5
2 6
3 5
1 7
5 7
7 8

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n, q = map(int, input().split())

# multiset< pair<int,int> > の代わりに、ソート済みリストで管理
# versions = [(version番号, 台数)]
versions = [(i, 1) for i in range(1, n + 1)]

for _ in range(q):
    x, y = map(int, input().split())
    res = 0
    new_versions = []

    # versions を小さい順に確認
    for v, cnt in versions:
        if v <= x:
            res += cnt  # アップグレード対象
        else:
            new_versions.append((v, cnt))  # そのまま残す

    # v ≤ x のものは削除され、y にまとめて追加される
    if res > 0:
        new_versions.append((y, res))

    # multiset と同じく、再度ソートしておく（O(M log M)）
    new_versions.sort()
    versions = new_versions

    print(res)

# note: このコード中の multiset を set にすると WA になります。何故でしょうか?
# 1️⃣ set は「重複要素」を許さない

# set には同じキー（例：(Y, 台数) の Y が同じ）が入れられません。
# しかしこの問題では、複数回のアップグレードによって同じバージョン番号に複数グループが存在することがあります。

# 2️⃣ 後で集約されるまでは、重複したままでよい

# multiset では (3,1), (3,2) のように「同じ version の異なるペア」が共存できます。
# 必要なときに v ≤ X の部分をまとめて削除・合算すればよく、途中で自動的にまとめる必要はありません。