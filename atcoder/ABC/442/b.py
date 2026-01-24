import io
import sys
from collections import Counter

# 下記に標準入力を記載
_InPUT = """\
10
2
1
3
1
3
1
1
3
2
2
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

q = int(input())

volume = 0
now_music = False
for i in range(q):
    query = int(input())
    if query == 1:
        volume += 1
    elif query == 2:
        if volume >= 1:
            volume -= 1
        # elif volume == 0:
            
    elif query == 3:
        if now_music == False:
            now_music = True
        elif now_music == True:
            now_music = False

    if volume >= 3 and now_music == True:
        print("Yes")
    else:
        print("No")



    
    