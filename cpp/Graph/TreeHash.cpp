#include <bits/stdc++.h>
using namespace std;
#define DEBUG 1
#define endl '\n'
#define FUCK if (DEBUG) cout << "fuck" << endl;
using ll = long long;
using pii = pair<int,int>;
using i128 = __int128_t;
using pll = pair<ll,ll>;
const ll N = 2000000;
const ll INF = 1e18;
const ll MOD = 1e9 + 7;

/*
定义 val[u] = 1 + f(h(val[v]))
v is son of u
*/

ll h(ll x) {
    return x * x * x * 1237123 + 19260817;
}
ll f(ll x) {
    ll cur = h(x & ((1 << 31) - 1)) + h(x >> 31);
    return cur;
}