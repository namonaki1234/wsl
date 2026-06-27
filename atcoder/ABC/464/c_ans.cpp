#include<bits/stdc++.h>
using namespace std;
using pi=pair<int,int>; // 古い色と新しい色のペアを扱いやすくするための準備

int main(){
  int n, m;
  cin >> n >> m;
  
  int kind = 0; // 現在の色の種類数
  vector<int> cnt(n+1); // cnt[k] は色kの鳥の数
  vector<vector<pi>> change(m+1); // change[d] は d日目に起きる変化のリスト
  
  // 【ステップ1＆2】初期状態の作成とスケジュール帳への記入
  for(int i=0; i<n; i++){
    int a, d, b;
    cin >> a >> d >> b;
    
    // 初期状態の記録
    if(cnt[a] == 0){ kind++; } // まだその色がなければ種類数を増やす
    cnt[a]++; // 色aの鳥の数を増やす
    
    // スケジュール帳に記入（d日目に a が減って b が増える）
    change[d].push_back({a, b}); 
  }
  
  // 【ステップ3】1日目からM日目までシミュレーション
  for(int i=1; i<=m; i++){
    // i日目の予定を順番に処理
    for(auto &nx : change[i]){
      // nx.first が古い色(減る色)、nx.second が新しい色(増える色)
      
      // 減る処理
      cnt[nx.first]--;
      if(cnt[nx.first] == 0){ kind--; }
      
      // 増える処理
      if(cnt[nx.second] == 0){ kind++; }
      cnt[nx.second]++;
    }
    
    // その日の種類数を出力
    cout << kind << "\n";
  }
  return 0;
}