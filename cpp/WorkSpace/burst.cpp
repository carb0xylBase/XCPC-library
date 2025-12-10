#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef pair<int,int> P;
#define rep(i,a,b) for(int i=(a);i<=(b);++i)
#define per(i,a,b) for(int i=(a);i>=(b);--i)
#define pb push_back
#define fi first
#define se second
#define SZ(x) (int)(x.size())
const int N=1e5+10,K=105;
const ll INF=1e10;
vector<P>e[N];
//0没有往上伸的一头 1有至少一条往上伸的一头 路径是可以共点的
int n,k,u,v,w,sz[N];
ll a[K],dp[N][K][2],tmp[K][2],ans,sum;
void upd(ll &x,ll y){
    x=max(x,y);
}
void dfs(int u,int fa){
    sz[u]=1;
    dp[u][0][0]=dp[u][1][1]=0;
    for(auto &x:e[u]){
        int v=x.fi,w=x.se;
        if(v==fa)continue;
        sum+=w;
        dfs(v,u);
        int up=min(sz[u],k);
        rep(i,0,up){
            rep(x,0,1){
                tmp[i][x]=dp[u][i][x];
            }
        }
        rep(i,0,up){
            for(int j=0;j<=sz[v] && i+j-1<=k;++j){
                rep(x,0,1){
                    rep(y,0,1){
                        if(tmp[i][x]<0 || dp[v][j][y]<0)continue;
                        upd(dp[u][i+j][x],tmp[i][x]+dp[v][j][y]);//没用
                        if(x==0){
                            if(y==1)upd(dp[u][i+j][1],tmp[i][x]+dp[v][j][y]+w);//用了 续在一头
                            upd(dp[u][i+j+1][1],tmp[i][x]+dp[v][j][y]+w);//多条路径
                        }
                        if(x==1){
                            if(y==1)upd(dp[u][i+j-1][0],tmp[i][x]+dp[v][j][y]+w);//用了 两头接在一起
                            upd(dp[u][i+j][0],tmp[i][x]+dp[v][j][y]+w);//用了 续在一头
                        }
                    }
                }
            }
        }
        sz[u]+=sz[v];
    }
    // rep(i,0,sz[u]){
    //     rep(j,0,1){
    //         if(dp[u][i][j]<0)continue;
    //         printf("u:%d i:%d j:%d dp:%d\n",u,i,j,dp[u][i][j]);
    //     }
    // }
}
int main(){
    freopen("input.txt", "r", stdin);
    int SSS = clock();
    scanf("%d%d",&n,&k);
    rep(i,1,n){
        rep(j,0,k){
            rep(x,0,1){
                dp[i][j][x]=-INF;
            }
        }
    }
    rep(i,2,n){
        scanf("%d%d%d",&u,&v,&w);
        e[u].pb(P(v,w));
        e[v].pb(P(u,w));
    }
    rep(i,1,k){
        scanf("%lld",&a[i]);
    }
    sort(a+1,a+k+1);
    rep(i,1,k){
        a[i]+=a[i-1];
    }
    dfs(1,0);
    //printf("sum:%d\n",sum);
    rep(i,0,k){
        rep(j,0,1){
            ans=max(ans,dp[1][i][j]-a[i]);
        }
    }
    printf("%lld\n",2*sum-ans);
    // cout << clock() - SSS << endl;
    return 0;
}