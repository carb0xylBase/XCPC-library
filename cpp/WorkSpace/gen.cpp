#include <bits/stdc++.h>
using namespace std;
using ll = long long;

int main() {
    cout << 1 << endl;
    srand(time(0));
    int n = 80;
    int m = rand() % n + 1;
    cout << n << ' ' << m << endl;
    for (int i = 1;i<m+1;i++) {
        cout << rand() % n + 1 << ' ';
    }
    cout << endl;
    for (int i = 2;i<=n;i++) {
        cout << i << ' ' << rand() % (i - 1) + 1 << endl;
    }
    return 0;
}
