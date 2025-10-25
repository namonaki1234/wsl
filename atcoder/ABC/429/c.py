import io
import sys
from collections import Counter
from itertools import combinations
import math

# 下記に標準入力を記載
_InPUT = """\
5
3 2 5 2 2
"""
sys.stdin = io.StringIO(_InPUT)
# 4
# 3 2 3 2
# 5
# 3 2 3 2 2
# ここからコードを記載

n = int(input())
a = list(map(int, input().split()))

# targets = []
# for i in range(n):
#     targets.append(a[i])
#     for j in range(i + 1, n):
#         if a[j] == targets[i]:
#             a[j] = -1

# a_counter = Counter(a)

# calc_counts = set()
# over_2_counts = set()
# for key in a_counter:
#     print(f"{key} : {a_counter[key]}")
#     if a_counter[key] >=2:
#         a_counter[key] = 2
#         over_2_counts.add(key)
#     calc_counts.add((key, a_counter[key]))

# print(calc_counts)
# print(over_2_counts)

index_counts = [0] * (max(a) + 1)

for i in a:
    index_counts[i] += 1

# print(index_counts)

# over_2count = []
# over_2_num = []
# count_1 = 0
# for i, count in enumerate(index_counts):
#     if count >= 2:
#         over_2count.append(count)
#         over_2_num.append(i)
#     elif count == 1:
#         count_1 += count

# print(over_2count)
ans = 0
# over_2count_combinations = list(combinations(over_2count, 2))
# over_2count_num_combinations = list(combinations(over_2_num, 2))
# over_2_num_count_combinations = [[x, y] for x, y in zip(combinations(over_2_num, 2), combinations(over_2count, 2))]
# print(over_2count_combinations)
# print(over_2count_num_combinations)
# print(over_2_num_count_combinations)
# if len(over_2count) == 0:
#     print(0)
#     exit()
# if count_1 >= 1 and len(over_2count) == 1:
#         ans += (math.comb((over_2count[0]), 2)) * count_1
# if count_1 >= 1 and len(over_2count) >= 2:
#     for i in range(len(over_2count)):
#         ans += (math.comb((over_2count[i]), 2)) * count_1
    
#     for i in range(len(over_2count_combinations)):
#         ans += over_2count_combinations[i][0] * over_2count_combinations[i][1]

# if len(over_2count) >= 2:
#     for i in range(len(over_2count_combinations)):
#             ans += over_2count_combinations[i][0] * over_2count_combinations[i][1]
#     ans += math.comb(len(over_2count), 2)
for i in range(len(index_counts)):
    if index_counts[i] >= 2:
        ans += math.comb(index_counts[i], 2) * (n - index_counts[i])


# for i in range(len(over_2count)):
#         ans += (math.comb((over_2count[i]), 2))
print(ans)