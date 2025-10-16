#include <bits/stdc++.h>
using namespace std;
int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
    int T = 1;
    int n = 5, m = 7;
    cout<<T<<"\n";
    for(int tc=0;tc<T;tc++){
        vector<pair<int,int>> edges;
        vector<int> p(n);
        iota(p.begin(),p.end(),1);
        shuffle(p.begin(),p.end(),rng);
        for(int i=1;i<n;i++) edges.emplace_back(p[i-1],p[i]);
        unordered_set<long long> seen;
        auto key=[&](int a,int b){ if(a>b) swap(a,b); return (long long)a<<32 | (unsigned long long)b; };
        for(auto &e:edges) seen.insert(key(e.first,e.second));
        uniform_int_distribution<int> dist(1,n);
        while((int)edges.size()<m){
            int u=dist(rng), v=dist(rng);
            if(u==v) continue;
            long long k=key(u,v);
            if(seen.count(k)) continue;
            seen.insert(k);
            edges.emplace_back(u,v);
        }
        cout<<n<<" "<<m<<"\n";
        for(auto &e:edges) cout<<e.first<<" "<<e.second<<"\n";
    }
}
