// g++ -D LOCAL c_ans.cpp
// ./a.out

#include <bits/stdc++.h>
using namespace std;

#define N (int)3e+5

int main(void){
// ==========================================
    // 提出時はここから下のブロックを削除・コメントアウトする
    // ローカル環境でのみ実行されるブロック
#ifdef LOCAL
    string test_input = R"(
3 7
1 1
1 3
1 3
2 1
2 2
1 2
2 1
)";
    stringstream ss(test_input);
    cin.rdbuf(ss.rdbuf()); // cin の読み込み先を test_input に変更
#endif
    // ==========================================

    // ここから下のコードは一切変更しなくてOK

	int a[N+1]={};
	int c[N+1]={};
	int n,q,t,x,mn=0;

	cin>>n>>q;
	for(int i=0;i<q;i++){
		cin>>t>>x;
		if(t==1){
			a[x]++;
			c[a[x]]++;
			if(c[a[x]]==n)mn=a[x];
		}
		if(t==2){
			if(x+mn>q)cout<<0<<endl;
			else cout<<c[x+mn]<<endl;
		}
	}

	return 0;
}
