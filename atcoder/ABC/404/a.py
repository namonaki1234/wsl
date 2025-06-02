from collections import defaultdict,deque
import io
import sys
import string


# 下記に標準入力を記載
_INPUT = """\
qazplwsxokmedcijnrfvuhbgt

"""


sys.stdin = io.StringIO(_INPUT)

S = input().strip()
alphabets = deque(string.ascii_lowercase)

for s in S:
    if s in alphabets:
        alphabets.remove(s)

print(alphabets[0])


