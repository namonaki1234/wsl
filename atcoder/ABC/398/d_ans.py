import io
import sys

# 標準入力を模擬するための入力データ
_INPUT = """\
6 -2 1
NNEEWS
"""
sys.stdin = io.StringIO(_INPUT)

# 入力の取得
n,t_x, t_y = map(int, input().split())  # 目標座標の取得
s = input().strip()  # 移動指示の文字列

# 座標をリストで管理（可変にするため）
target = [t_x, t_y]  # 目標座標
current = [0, 0]     # 現在の座標

# 訪問済み座標の集合（タプルで管理）
visited = set()
visited.add((current[0], current[1]))  # 初期位置を訪問済みに追加

result = []  # 結果を格納するリスト

# 移動指示に従って処理
for move in s:
    # 移動処理
    if move == 'N':
        target[0] += 1
        current[0] += 1
    elif move == 'W':
        target[1] += 1
        current[1] += 1
    elif move == 'S':
        target[0] -= 1
        current[0] -= 1
    elif move == 'E':
        target[1] -= 1
        current[1] -= 1
    
    # 現在の座標を訪問済みとして記録(setだから逐一全部ぶちこんでも重複は無視されて、新たな煙としてカウントされる)
    visited.add((current[0], current[1]))
    
    # 現在の目標座標が訪問済みかチェック
    if (target[0], target[1]) in visited:
        result.append('1')  # 訪問済みなら'1'を追加
    else:
        result.append('0')  # 訪問していなければ'0'を追加

# 結果を出力
print(''.join(result))
