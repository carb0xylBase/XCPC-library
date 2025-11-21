#define DEBUG 1
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

ll L[N], R[N], b[N];
ll tpL = -1, tpR = -1, tpb = -1;

void solve() {
    int n; cin >> n;
    vector<ll> a(n + 1, 0);
    for (int i = 1;i<n+1;i++) cin >> a[i];
    int m; cin >> m;
    while (m --) {
      ll l, r, x; cin >> l >> r >> x;
      // vector<ll> b;
      tpb = -1;
      for (int i = l;i<=r;i++) b[++tpb] = a[i];
      ll cnt = 0, sum = 0;
      ll ans = 0;
      while (tpb != -1) {
        bool ok = 0;
        if (tpb == 0) ok = 1;
        nth_element(b, b + (tpb + 1) / 2, b + tpb + 1);
        ll val = b[(tpb + 1) / 2];
        ll M = 0;
        tpL = tpR = -1;
        for (int i = 0;i<=tpb;i++) {
          ll v = b[i];
          if (v < val) {
            L[++tpL] = v;
            // L.push_back(v);
          } else if (v > val) {
            // R.push_back(v);
            R[++tpR] = v;
          } else {
            M ++;
          }
        }
        while (M --) {
          if (tpL <= tpR) L[++tpL] = val;
          else R[++tpR] = val;
        }
        ll res = cnt * val - sum;
        for (int i = 0;i<=tpL;i++) {
          ll v = L[i];
          res += val - v;
        }
        if (res <= x) {
          ans = val;
          for (int i = 0;i<=tpL;i++) {
            ll v = L[i];
            cnt ++;
            sum += v;
          }
          for (int i = 0;i<=tpR;i++) {
            b[i] = R[i];
          }
          tpb = tpR;
        } else {
          // b.swap(L);
          for (int i = 0;i<=tpL;i++) {
            b[i] = L[i];
          }
          tpb = tpL;
        }
        if (ok) break;
      }

      cnt = 0;
      for (int i = l;i<=r;i++) {
        if (a[i] <= ans) {
          cnt ++;
          x -= ans - a[i];
          a[i] = ans;
        }
      }
      ll yu = x - x / cnt * cnt;
      for (int i = l;i<=r;i++) {
        if (a[i] <= ans) {
          a[i] += x / cnt;
          if (yu) {
            yu --;
            a[i] ++;
          }
        }
      }
    }

    for (int i = 1;i<=n;i++) {
      cout << a[i] << ' ';
    }
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
