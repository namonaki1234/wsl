import io
import sys
from collections import Counter,deque

# 下記に標準入力を記載
_InPUT = """\
3
3 1
7 2 3
1 3 2
3 2
7 2 3
1 3 2
2 1
2 1
1 2
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
#卵は 在庫にある卵の中で仕入れた順に使用される→queueを使う
#n日営業
#i日目朝にA_i個の卵を仕入れる。昼にB_iの卵を使用。夜にD日以上経過した卵を処分する
#n日後の夜時点で何個の卵が残っているのか
"""
この問題で一番大事なのは、
キューに「卵の個数」を入れるのではなく，「その卵を仕入れた日」を入れる

今のコードだとq.append(a[i])
になっていますが，これは「i日目に a[i] 個仕入れた」という情報を
a[i] という1個の数字で入れてしまっている

でも本当に知りたいのは，
この卵は何日目に仕入れたか
いま一番古い卵は何日目のものか
"""
# 卵は仕入れた順に使う → queue
# キューには「その卵を仕入れた日」を入れる

t = int(input())

for _ in range(t):
    n,d = map(int,input().split())
    a = [int(x) for x in input().split()]
    b = [int(x) for x in input().split()]

    q = deque()

    for i in range(n):
        # 朝：i日目に a[i] 個仕入れる
        for _ in range(a[i]):
            q.append(i)
        # 昼：古い順に b[i] 個使う
        for _ in range(b[i]):
            q.popleft()
        
         # 夜：D日以上経った卵を捨てる
         #ヒープの時と同様先頭だけ見れば良い（卵は古い順に入っているから）
        while q and q[0] == i - d:
            q.popleft()

    print(len(q))            
