import io
import sys
from collections import defaultdict,deque,Counter

# 下記に標準入力を記載
_InPUT = """\
17
AAABABABBBABABBABABABABBAAABABABBA
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# sys.setrecursionlimit(1000000)

n = int(input())
s = input()

len_s = len(s)
kinds = sorted(set(s))

cnt = Counter(s)
a, b = kinds[0], kinds[1]
ca, cb = cnt[a], cnt[b]

# 交互並びの必要条件（個数差が0）
# if abs(ca - cb) != 0:
#     print(-1)
#     exit()

# a の現在位置配列
pos_a = []
for i, ch in enumerate(s):
    if ch == a:
        pos_a.append(i)

# 偶数インデックスのスロット数
need_even = (len_s + 1) // 2  
 # 奇数インデックスのスロット数
need_odd  = len_s // 2       

INF = 10**18
ans = INF

# パターン1: a を偶数番目（0,2,4,...) に配置できる個数条件ならコスト計算
if ca == need_even and cb == need_odd:
    # targets = [0,2,4,...]
    cost = 0
    t = 0
    for p in pos_a:
        cost += abs(p - t)
        t += 2
    ans = min(ans, cost)

# パターン2: a を奇数番目（1,3,5,...) に配置できる個数条件ならコスト計算
if ca == need_odd and cb == need_even:
    # targets = [1,3,5,...]
    cost = 0
    t = 1
    for p in pos_a:
        cost += abs(p - t)
        t += 2
    ans = min(ans, cost)

print(ans)


# def ok(t: str) -> bool:
#     # どの隣り合う位置でも同じ文字が隣接しないか
#     for i in range(len(t) - 1):
#         if t[i] == t[i+1]:
#             return False
#     return True

# # すでに条件を満たしていれば0
# if ok(s):
#     print(0)
#     exit()

# visited = set()
# visited.add(s)

# dist = {s: 0}

# Q = deque([s])

# while Q:
#     pos = Q.popleft()
#     d = dist[pos]

#     #隣接スワップを全通り試す
#     for i in range(2*n-1):
#         s_list = list(pos)
#         s_list[i], s_list[i+1] = s_list[i+1], s_list[i]
#         new_s = ''.join(s_list)

#         if new_s in visited:
#             continue
    
#         visited.add(new_s)
#         dist[new_s] = d + 1
#         if ok(new_s):
#             print(dist[new_s])
#             exit()
#         Q.append(new_s)

