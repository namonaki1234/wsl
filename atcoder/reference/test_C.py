import io
import sys
from collections import deque

# 標準入力を設定
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

# クエリの数を取得
Q = int(input())
queries = [list(map(int, input().split())) for _ in range(Q)]

# ヘビの待ち行列とオフセットの初期化
queue = deque()
offset = 0
results = []
var = []
l = deque()

# クエリを順次処理
for query in queries:
    query_type = query[0]
    if query_type != 2: 
        variable = query[1]
        var.append(variable)
    
    if query_type == 1:  # タイプ1のクエリ
        length = query[1]
        # 新しいヘビの頭の位置を計算
        if not queue:  # 列が空の場合
            queue.append(length)  # 最初のヘビの頭の位置
            l.append(length)
        else:  # 列にヘビがいる場合
            last_head = queue[-1]  # 最後尾のヘビの頭の位置
            queue.append(last_head + length)  # 新しいヘビの頭の位置を追加
            l.append(length)

    elif query_type == 2:  # タイプ2のクエリ
            removed_length = l.popleft()  # 先頭のヘビを列から抜けさせる
            removed_queue = queue.popleft()
            offset += removed_length  # オフセットを更新

    elif query_type == 3:  # タイプ3のクエリ
        k = variable - 1
        head_position = queue[k] - offset - l[k] # オフセットを考慮して頭の位置を計算
        print(head_position)
        # results.append(head_position)  # 結果をリストに追加


