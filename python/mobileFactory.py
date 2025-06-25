# # 日をまたぐ時間範囲をチェックする関数
# def is_in_range(target, start, end):
#     if start == end:
#         return True  # 同じなら含むとする
#     if start < end:
#         return start <= target < end
#     else:
#         # 日をまたぐ場合（例：22〜5）
#         return target >= start or target < end

# # テスト例
# print(is_in_range(3, 22, 5))  # True
# print(is_in_range(6, 22, 5))  # False
# print(is_in_range(22, 22, 5)) # True
# print(is_in_range(4, 4, 4))   # True

# import io
# import sys

# # 下記に標準入力を記載
# _InPUT = """\
# 3
# 22
# 5
# """
# sys.stdin = io.StringIO(_InPUT)
# # ここからコードを記載

print("ある時刻（0~23）を入力してください")
target = int(input())
print("開始時刻（0~23）を入力してください")
start = int(input())
print("終了時刻（0~23）を入力してください")
end = int(input())

#例外は一番最初のif文で処理する
if start == end:
    print("同じ時刻の場合は範囲に含むと判断する。")
else:
    if start < end:
        is_in_range = start <= target < end
    else:
        # 日をまたぐ場合（例：22〜5）
        is_in_range = target >= start or target < end

    if is_in_range:
        print(f"{target}時は{start}時から{end}時の範囲に含まれる。")
    else:
        print(f"{target}時は{start}時から{end}時の範囲に含まれない。")
