#include <bits/stdc++.h>
using namespace std;
using ll = long long;
int main(){
    mt19937_64 rng((unsigned)chrono::high_resolution_clock::now().time_since_epoch().count());
    uniform_int_distribution<int> coord(-10,10);
    uniform_int_distribution<int> rad(1,20);
    int N=2;
    struct C{int x,y,r;};
    vector<C> a;
    // make first big circle
    C first;
    first.r = uniform_int_distribution<int>(10,20)(rng);
    first.x = coord(rng);
    first.y = coord(rng);
    a.push_back(first);
    for(int i=1;i<N;i++){
        int tries=0;
        while(tries<2000){
            int r = uniform_int_distribution<int>(1, min(19, a[0].r-1))(rng);
            int x = uniform_int_distribution<int>(-10,10)(rng);
            int y = uniform_int_distribution<int>(-10,10)(rng);
            bool ok=true;
            for(auto &c:a){
                long long dx=c.x-x, dy=c.y-y;
                long long d2=dx*dx+dy*dy;
                long long rr=c.r - r;
                if(rr>=0){
                    if(d2>rr*rr) ok=false;
                }else{
                    // new circle can't contain existing one, avoid intersection: require separation
                    long long sum=c.r + r;
                    if(d2<=(long long)sum*sum) ok=false;
                }
                if(!ok) break;
            }
            if(ok){
                a.push_back({x,y,r});
                break;
            }
            tries++;
        }
        if((int)a.size()==i){ // fallback: place strictly inside first
            int r = max(1, a[0].r - i - 1);
            int x = a[0].x;
            int y = a[0].y;
            a.push_back({x,y,r});
        }
    }
    cout<<a.size()<<"\n";
    for(auto &c:a) cout<<c.x<<" "<<c.y<<" "<<c.r<<"\n";
    return 0;
}
