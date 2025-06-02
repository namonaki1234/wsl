import io
import sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')
# 下記に標準入力を記載
_INPUT = """\
5 2 4 1 2 11 3 
10
"""

sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載
# 巻き尺おじさんのお仕事
# coding: utf-8
# Your code here!
li = list(map(int, input().split()))
X = int(input())
# print(li, X)
N = len(li)

right, ans, measure = 0, 0, 0

for left in range(N):
    print("left詰める","l:r",left,right,"m",measure)
    while right < N and measure+li[right] <= X:
        measure += li[right]
        right += 1
        print("-rightのびる","l:r",left,right,"m",measure)
    ans = max(ans, right-left)
    if left==right:
        right += 1
        continue
    measure -= li[left]
print(ans)
        