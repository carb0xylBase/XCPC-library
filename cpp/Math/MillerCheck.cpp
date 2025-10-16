#include <bits/stdc++.h>
using namespace std;
using ull = unsigned long long;
using u128 = __uint128_t;

ull modmul(ull a, ull b, ull m){ return (u128)a*b % m; }
ull modpow(ull a, ull e, ull m){ u128 r=1, x=a%m; while(e){ if(e&1) r=(r*x)%m; x=(x*x)%m; e>>=1; } return (ull)r; }
bool miller_check(ull n, ull a){
    if(a % n == 0) return true;
    ull d = n-1; int s = 0;
    while((d&1)==0){ d >>= 1; ++s; }
    ull x = modpow(a, d, n);
    if(x==1 || x==n-1) return true;
    for(int i=1;i<s;i++){
        x = modmul(x, x, n);
        if(x==n-1) return true;
    }
    return false;
}
// 时间复杂度: O(k * log^3 n)（k为基础集合大小，常数）
// 用法: if(isPrime(n)) // n <= 1e18
bool isPrime(ull n){
    if(n < 2) return false;
    for(ull p : {2ull,3ull,5ull,7ull,11ull,13ull,17ull,19ull,23ull,29ull,31ull,37ull}){
        if(n == p) return true;
        if(n % p == 0) return false;
    }
    ull bases[] = {2ull,325ull,9375ull,28178ull,450775ull,9780504ull,1795265022ull};
    for(ull a : bases) if(!miller_check(n, a)) return false;
    return true;
}
