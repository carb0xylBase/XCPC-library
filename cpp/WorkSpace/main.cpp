#include <bits/stdc++.h>
using namespace std;
#define DEBUG 1
#define endl '\n'
#define FUCK if (DEBUG) cout << "fuck" << endl;
using ll = long long;
using pii = pair<int,int>;
using pll = pair<ll,ll>;
using db = double;
using pdd = pair<double, double>;
using i128 = __int128_t;
const ll N = 2000000;
const ll INF = 1e18;
const ll MOD = 1e9 + 7;

void solve() {
    
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
    cin >> t;

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