import io
import sys
from collections import Counter

# 下記に標準入力を記載
_InPUT = """\
beginner
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# lはリストの長さ（選べるジュースの数）、xはジュースのid識別番号

s = input()

s_counter = Counter(s)
# print(s_counter)
# print(s_counter.values())
# print(s_counter.most_common(2))

max_cnt = max(s_counter.values())
# print(max_cnt)

s_max_cnt = [k for k,v in s_counter.items() if v == max_cnt]

# print(s_max_cnt)
# s.remove("i")
s_remove = [x for x in s if x not in s_max_cnt]

print("".join(s_remove))