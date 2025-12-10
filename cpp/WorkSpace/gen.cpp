#define DEBUG 0
#define FUCK cout << "fuck" << endl;
#if DEBUG
    #include "all.hpp"
#else
    #include <bits/stdc++.h>
#endif

using namespace std;
using ll = long long;
using pii = pair<int,int>;
using pll = pair<ll,ll>;
using db = long double;
using pdd = pair<db, db>;
using i128 = __int128_t;

const ll N = 2000000;
const ll INF = 5e18;
const ll MOD = 1e9 + 7;

void solve() {
    srand(time(0));
    // freopen("input.txt", "w", stdout);
    int n = 6;
    int K = 2;
    cout << n << ' ' << K << endl;
    for (int i = 2;i<=n;i++) {
        cout << i << ' ' << rand() % (i - 1) + 1 << ' ' << rand() % 10 + 1 << endl;
    }
    // for (int i = 2;i<=n;i++) {
    //     cout << i - 1 << ' ' << i << ' ' << 100 << endl;
    // }
    for (int i = 1;i<=K;i++) {
        cout << rand() % 10 + 1 << ' ';
    }
    cout << endl;
    return;
}

signed main() {
#if DEBUG
    freopen("input.txt", "r", stdin);
    auto start_time = chrono::steady_clock::now();
#else
    ios::sync_with_stdio(false);
#endif
    cin.tie(nullptr);

    int t = 1;
    // cin >> t;

    while (t--) {
        solve();
    }

#if DEBUG
    auto end_time = chrono::steady_clock::now();
    auto diff = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
    cerr << "Time: " << diff.count() << " ms" << endl;

#endif

    return 0;
}

/*
4 4
2 1 9
3 1 5
4 1 4
1 8 10 3 

*/