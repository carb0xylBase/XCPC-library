#include <bits/stdc++.h>
using namespace std;
#define DEBUG 1
using ll = long long;
using pii = pair<int,int>;
using pll = pair<ll,ll>;
const ll N = 4000000;
const ll INF = 5e18;
const ll MOD = 1e9 + 7;

namespace LagInt {
    ll qpow(ll base,ll k,ll mod) {
        ll res = 1;
        if (k < 0) {
            return qpow(qpow(base,-k,mod),mod-2,mod);
        }
        while (k) {
            if (k & 1) {
                res *= base; res %= mod;
            }
            k >>= 1;
            base *= base; base %= mod;
        }
        return res;
    }

    ll fac[N],facInv[N];

    void init(int n) {
        fac[0] = 1; facInv[0] = 1;
        for (int i = 1;i<=n;i++) {
            fac[i] = fac[i - 1] * i % MOD;
        }
        
        vector<ll> pre(n + 1,1),suc(n + 1,1);
        for (int i = 1;i<=n;i++) {
            pre[i] = pre[i - 1] * fac[i] % MOD;
        }
        suc[n] = fac[n];
        for (int i = n - 1;i>=0;i--) {
            suc[i] = suc[i + 1] * fac[i] % MOD;
        }

        ll S = qpow(pre[n],MOD - 2,MOD);

        for (int i = 1;i<=n;i++) {
            facInv[i] = S * pre[i - 1] % MOD * suc[i + 1] % MOD;
        }
        return;
    }

    // 当 x[i] = i 时,用 A 模式,复杂度 nlogn
    // 记得先调用预处理
    ll LagrangeInterpolationA(ll x,vector<ll>& a) {
        x %= MOD;
        if (x < 0) return 0;
        if (a.size() > x) {
            return a[x];
        }
        ll t = 1;
        for (int i = 1;i<=a.size() - 1;i++) {
            t *= (x - i); t %= MOD;
        }

        ll ans = 0;
        int n = int(a.size()) - 1;

        for (int i = 1;i<a.size();i++) {
            ll y = a[i];
            ll res = t;
            res *= y; res %= MOD;
            res *= qpow((x - i) % MOD,MOD-2,MOD); res %= MOD;
            res *= facInv[i - 1]; res %= MOD;
            res *= facInv[n - i]; res %= MOD;
            ll sign = ((n - i) & 1) ? (MOD - 1) : 1;
            res = res * sign % MOD;
            ans += res; ans %= MOD;
        }
        return ans;
    }
};

void solve() {
    ll n,k;
    cin >> n >> k;
    LagInt::init(2e6);

    vector<ll> a(k + 3,0);
    for (int i = 1;i<k+3;i++) {
        a[i] = a[i - 1] + LagInt::qpow(i,k,MOD);
        a[i] %= MOD;
    }

    cout << LagInt::LagrangeInterpolationA(n,a) << endl;
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
