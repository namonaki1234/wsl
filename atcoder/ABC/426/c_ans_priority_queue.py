import io
import sys
import heapq

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

# priority_queue< pair<int,int>, vector<>, greater<> >
# → Python では heapq（最小ヒープ）で実現
pq = []
for i in range(1, n + 1):
    heapq.heappush(pq, (i, 1))  # (バージョン番号, 台数)

for _ in range(q):
    x, y = map(int, input().split())
    res = 0

    # バージョン番号が x 以下のものを全部取り出す
    while pq and pq[0][0] <= x:
        v, cnt = heapq.heappop(pq)
        res += cnt  # 今回アップグレードされた PC 台数に加算

    # アップグレードされたものを y バージョンとして戻す
    if res > 0:
        heapq.heappush(pq, (y, res))

    print(res)
