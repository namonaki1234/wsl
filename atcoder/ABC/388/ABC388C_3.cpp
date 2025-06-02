#include <iostream>
#include <vector>
using namespace std;

int main() {
    int N;
    cin >> N;

    vector<int> A(N);
    for (int i = 0; i < N; ++i) {
        cin >> A[i];
    }

    int kind_count = 0;


    if (A[N - 1] / 2 < A[0]) {
        cout << 0 << endl;
    } else {
        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < N; ++i) {
                if (A[N - 1 - j] / 2 >= A[i]) {
                    kind_count++;
                }
            }
        }
        cout << kind_count << endl;
    }

    return 0;
}
