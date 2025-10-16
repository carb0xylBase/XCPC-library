#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
#define DEBUG 1
const ll N = 2000000;

ll qpow(ll base,ll k,ll mod) {
    if (base == 0) return 0;
    ll res = 1;
    base %= mod; base = (base + mod) % mod;
    k %= (mod - 1); k = (k + mod - 1) % (mod - 1);
    while (k) {
        if (k & 1) {    
            res *= base; res %= mod;
        }
        k >>= 1;
        base *= base; base %= mod;
    }
    return res;
}

ll phi(ll x) {
    ll y = x;
    vector<ll> q;
    for (int i = 2;i*i<=y;i++) {
        if (y % i == 0) {
            q.push_back(i);
        }

        while (y % i == 0) {
            y /= i;
        }
    }

    if (y > 1) q.push_back(y);

    for (auto v : q) {
        x /= v;
        x *= (v - 1);
    }

    return x;
}

class Primes{
public:
    ll notPrime[N];
    ll phi[N],mu[N];
    vector<ll> primes;
    void sieve(int maxn){
        phi[1] = 1;
        mu[1] = 1;
        for (ll i = 2;i<maxn + 1;i++){
            if (!notPrime[i]){
                primes.push_back(i);
                phi[i] = i - 1;
                mu[i] = -1;
            }

            for (auto p : primes) {
                if (i * p > maxn) break;
                notPrime[i * p] = 1;
                if (i % p == 0) {
                    phi[i * p] = phi[i] * p;
                    mu[i * p] = 0;
                    break;
                }
                phi[i * p] = phi[i] * phi[p];
                mu[i * p] = -mu[i];
            }
        }
        return;
    }
};

// 找奇质数 x 的原根,一定是奇素数!!!
// 时间复杂度不详
ll getYuanGen(ll x) {
    Primes solver;
    solver.sieve(x);
    ll y = x - 1;
    vector<ll> q;
    for (int i = 0;i<solver.primes.size() && y > 1;i++) {
        if (y % solver.primes[i] == 0) {
            q.push_back(solver.primes[i]);
        }
        while (y % solver.primes[i] == 0) {
            y /= solver.primes[i];
        }
    }

    for (ll i = 1;i<x;i++) {
        bool ok = 1;
        for (auto v : q) {
            if (qpow(i,(x-1)/v,x) == 1) {
                ok = 0;
                break;
            } else {

            }
        }
        if (ok) return i;
    }

    return 0;
}