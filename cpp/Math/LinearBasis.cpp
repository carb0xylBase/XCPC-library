#include <bits/stdc++.h>
using ll = long long;
const ll N = 2000000;
const ll INF = 5e18;
const ll MOD = 1e9 + 7;
using namespace std;

struct LinearBasis {
    ll basis[63];
    bool flag0; // 是否可以表示 0

    void insert(ll x) {
        if (!x) return;
        for (int i = 60;i>=0;i--) {
            if ((1ll<<i)&x) {
                if (basis[i]) {
                    x ^= basis[i];
                } else {
                    basis[i] = x;
                    break;
                }
            }
        }
        if (x) flag0 = 1;
        return;
    }

    void preprocess() {
        for (int i = 60;i>=0;i--) {
            if (!basis[i]) continue;
            for (int j = i-1;j>=0;j--) {
                if (basis[j] && ((1ll<<j)&basis[i])) {
                    basis[i] ^= basis[j];
                }
            }
        }
        return;
    }

    void init() {
        flag0 = 0;
        memset(basis,0,sizeof basis);
        return;
    }

    ll get_mn() {
        for (int i = 0;i<=60;i++) {
            if (basis[i]) {
                return basis[i];
            }
        }
    }

    ll get_mx() {
        ll res = 0;
        for (int i = 60;i>=0;i--) {
            if (basis[i] && !((1ll<<i)&res)){
                res ^= basis[i];
            }
        }
        return res;
    } 
} ;