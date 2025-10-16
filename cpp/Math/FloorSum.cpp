#include <bits/stdc++.h>
using namespace std;
using i128 = __int128_t;
using ll = long long;

// floor((a * i + b) / m)
// i 从 0 到 n - 1

i128 floor_sum_i128(i128 n, i128 m, i128 a, i128 b) {
    i128 ans = 0;
    if (n <= 0) return 0;
    while (true) {
        if (a >= m) {
            ans += (a / m) * (n * (n - 1) / 2);
            a %= m;
        }
        if (b >= m) {
            ans += (b / m) * n;
            b %= m;
        }
        i128 y = a * n + b;
        if (y < m) break;
        n = y / m;
        b = y % m;
        i128 tmp = m;
        m = a;
        a = tmp;
    }
    return ans;
}

ll modmul_i128(i128 x, ll mod) {
    ll r = (ll)(x % mod);
    if (r < 0) r += mod;
    return r;
}

ll qpow(ll a, ll e, ll mod) {
    ll r = 1 % mod;
    while (e) {
        if (e & 1) r = (i128)r * a % mod;
        a = (i128)a * a % mod;
        e >>= 1;
    }
    return r;
}

ll floor_sum_mod(i128 n, i128 m, i128 a, i128 b, ll mod) {
    if (n <= 0) return 0;
    ll ans = 0;
    ll inv2 = (mod % 2 == 1) ? qpow((mod + 1) / 2, 1, mod) : -1;
    while (true) {
        if (a >= m) {
            ll k = (ll)(a / m % mod);
            i128 t = n * (n - 1) / 2;
            ll tmod;
            if (inv2 != -1) {
                ll nm = (ll)(n % mod);
                ll nm1 = (ll)((n - 1) % mod);
                tmod = (i128)nm * nm1 % mod * inv2 % mod;
            } else {
                tmod = (ll)(t % mod);
            }
            ans = (ans + (i128)k * tmod) % mod;
            a %= m;
        }
        if (b >= m) {
            ll k = (ll)(b / m % mod);
            ll nmod = (ll)(n % mod);
            ans = (ans + (i128)k * nmod) % mod;
            b %= m;
        }
        i128 y = a * n + b;
        if (y < m) break;
        i128 n2 = y / m;
        i128 b2 = y % m;
        ll n2mod = (ll)(n2 % mod);
        ans = (ans + n2mod * (ll)(n % mod)) % mod;
        n = n2;
        b = b2;
        i128 tmp = m;
        m = a;
        a = tmp;
    }
    if (ans < 0) ans += mod;
    return ans;
}
