import io
import sys

# 下記に標準入力を記載
_InPUT = """\
35 3
thequickbrownfoxjumpsoverthelazydog

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n,k = map(int, input().split())
s = input().strip()

max_count = 0
count_and_target = set()
for i in range(n-k+1):
    target = s[i:i+k]
    count = 0
    for j in range(n-k+1):
        if target == s[j:j+k]:
            count += 1
    max_count = max(count, max_count)
    count_and_target.add((count, target))

print(max_count)
count_and_target = list(count_and_target)
count_and_target.sort()
ans_strings = []
for c, t in count_and_target:
    if c == max_count:
        ans_strings.append(t)
print(' '.join(ans_strings))
