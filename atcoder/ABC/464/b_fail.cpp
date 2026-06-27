#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int main() {
    int h,w,count=0;
    cin >> h >> w;
    vector<string> c(h);

    for (int i=0;i<h;i++){
        cin >> c[i];
    }

    vector<string> c_copy = c;

    for (int i=h-1;i>=0;i--){
        if (c[i].find('#')!= string::npos){
                break;
            }
        count = 0;
        for (int j=0;j<w;j++){
            if (c[i][j]=='.'){
                count += 1;
            }
        }
        if (count == w){
            c_copy.erase(c_copy.begin() + i);
        }
        cout << count ;
    }


    for (int i=w-1;i>=0;i--){
        if (c[i].find('.')!= string::npos){
                break;
            }
        count = 0;
        for (int j=0;j<h;j++){
            if (c[j][i]=='.'){
                count += 1;
            }
        }
        if (count == h){
            for (int k=0;k<c_copy.size();k++){
                c_copy[k].erase(i,1);
            }
        }
        cout << count ;
    }

    
    for (int i=0;i<c_copy.size();i++){
        cout << c_copy[i] << endl;
    }

    // for (int i=0;i<c.size();i++){
    //     cout << c[i] << endl;
    // }
    
    
    return 0;
}