#include <bits/stdc++.h>
using namespace std;
#define DEBUG 1
using ll = long long;
using pii = pair<int,int>;
using i128 = __int128_t;
using pll = pair<ll,ll>;
const ll N = 2000000;
const ll INF = 5e18;
const ll MOD = 1e9 + 7;

template<ll M>
struct AutoMod {
    ll v;
    AutoMod(ll x=0){ v = x%M; if(v<0) v+=M; }
    static ll qpow(ll a,ll b){
        ll r=1, m=M;
        if(b<0) return qpow(qpow(a,-b),M-2);
        while(b){ if(b&1) r=r*a%m; a=a*a%m; b>>=1; }
        return r;
    }
    AutoMod pow(ll k)const{ return qpow(v,k); }
    AutoMod inv()const{ return qpow(v,M-2); }
    AutoMod operator+(AutoMod o) const{ return AutoMod(v+o.v); }
    AutoMod operator-(AutoMod o) const{ return AutoMod(v-o.v); }
    AutoMod operator*(AutoMod o) const{ return AutoMod(v*o.v); }
    AutoMod operator/(AutoMod o) const{ return *this * o.inv(); }
    bool operator==(AutoMod o) const{ return v==o.v; }
};
