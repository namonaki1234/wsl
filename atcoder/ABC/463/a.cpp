// g++ -D LOCAL c_ans.cpp
// ./a.out

#include <iostream>
#include <string>
using namespace std;

int main() {
    int x,y;
    cin >> x >> y;

    if ((x * 9)== (y * 16) ) {
        cout << "Yes";
    }
    else {
        cout << "No";
    }

    return 0;
}