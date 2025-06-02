from collections import defaultdict
import io
import sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')

# 下記に標準入力を記載
_INPUT = """\
40.0

"""


sys.stdin = io.StringIO(_INPUT)

X = float(input())

if X >= 38.0:
    print(1)
elif X >= 37.5:
    print(2)
else:
    print(3)


