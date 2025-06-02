import io
import sys
from collections import deque

# 下記に標準入力を記載
_INPUT = """\
7
1 5
1 7
3 2
1 3
1 4
2
3 3


"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

Q = int(input())
x = [list(map(int, input().split())) for l in range(Q)]
# クエリの数
queue = deque()  # ヘビの待ち行列
offset = 0  # オフセット
results = []  # 出力結果を格納するリスト
l = []  # ヘビの長さ



# クエリを順次処理
# for i in range(1, Q + 1):
    # query = list(map(int, x[i].split()))  # クエリを整数に変換
for i in range( Q ):
    if x[i][0] == 1:  # タイプ1のクエリ
        l.append(x[i][1]) # ヘビの長さ
    
        if not queue:  # 列が空の場合
            queue.append(l)  # 最初のヘビの尻尾の位置(次の頭の位置)を0に設定
        else:  # 列にヘビがいる場合
            last_head = queue[-1]  # 最後尾のヘビの尻尾の位置
            queue.append(last_head + l )  # 新しいヘビの頭の位置を追加
    elif x[i][0] == 2:  # タイプ2のクエリ
        m = queue.popleft()  # 先頭のヘビを列から抜けさせる
        offset += m  # オフセットを更新
    elif x[i][0] == 3:  # タイプ3のクエリ
        k = x[i][1] - 1  # 0-indexedに変換
        head_position = queue[k] - offset  # オフセットを考慮して頭の位置を計算
        results.append(head_position)  # 結果をリストに追加

# 結果を出力
print('\n'.join(map(str, results)))


