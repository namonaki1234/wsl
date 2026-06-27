// g++ -D LOCAL c_ans.cpp
// ./a.out

#include <iostream>
#include <string>
using namespace std;

int main() {
    string s;
    cin >> s;
    int count_e = 0;
    
    for (int i=0;i<s.size();i++){
        if (s[i] == 'E'){
            count_e += 1;
        }
    }

    // cout << count_e;
    if (count_e >= ((s.size()+1)/2) ) {
        cout << "East";
    }
    else {
        cout << "West";
    }

    return 0;
}