from collections import defaultdict
import io
import pprint
import sys
import math

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
d = defaultdict(int)
test_dict = dict(a=1, b=2, c=3)
print(list(test_dict.keys()))
print(test_dict.values())
print(test_dict.items())
print(test_dict.get('a'))
test_dict.setdefault('d', 4)
test_dict['a'] = 100
print(test_dict.items())

f = math.pi
print(round(f, 2))

for i,j in zip(range(10), range(10)):
    d[(i, j)] += 1

l = [ i for i in range(10)]
m = [ i+3 for i in range(10)]
a_list = [ i for i in range(10)]
a_tuple = tuple(a_list)
d_2 = dict(zip(l, m))
print(d_2)
n =[1, 2, 3, 4, 4,8,7,2]
n.insert(2, 1000)
n.extend([100, 200])
n.remove(1000)
n.pop(0)
# n.clear()
print(n)
set_n = set(n)
print(set_n)

for i,j in enumerate(l,start=1):
    # startを入れることでlの最初のindexを0から2（任意の数）に変更できる
    print(i, j,end=' , ')

for i in zip(l, m):
    print(i, end=' , ')

# zip_list = list(zip(l, m))

# pprint.pprint(zip_list,indent = 5, width=20)

# df = pd.DataFrame(zip_list, columns=['l', 'm'])
# print(N, M)
# print(d)
# print(id(d))