#include <iostream>
#include <cstring>
#include <ctime>
#include <algorithm>
#include <cmath>
#include <vector>
#include <queue>
#include <set>
#include <map>
using namespace std;
typedef long long int ll;
const int N=1e6+9,INF=1e9;
const double eps=1e-5;
typedef pair <int,int> PII;
inline int read()
{
    int x=0,f=1;char c=getchar();
    while(c<'0' || c>'9') {if(c=='-') f=-1;c=getchar();}
    while(c>='0' && c<='9') {x=x*10+c-48,c=getchar();}
    return x*f;
}
inline ll readll()
{
    ll x=0,f=1;char c=getchar();
    while(c<'0' || c>'9') {if(c=='-') f=-1;c=getchar();}
    while(c>='0' && c<='9') {x=x*10+c-48,c=getchar();}
    return x*f;
}
struct node{
    int to,nxt;
} edge[N];
int n,m,tot,head[N],siz[N],ans;
bool vis[N];
void addedge(int u,int v)
{
    edge[++tot].to=v,edge[tot].nxt=head[u],head[u]=tot;
}
void dfs(int u,int fa)
{
    for(int i=head[u];i;i=edge[i].nxt)
    {
        int v=edge[i].to;
        if(v==fa) continue;
        dfs(v,u);
        siz[u]+=siz[v];
    }
    if(vis[u])
    {
        ans+=(siz[u]+1)/2;
        if(siz[u]&1 && fa!=0) vis[fa]=true; 
        siz[u]=0;
    }
    else siz[u]++;
}
int main()
{
    // #define FILEIO
    #ifdef FILEIO
        freopen("in.in","r",stdin);
        freopen("out.out","w",stdout);
    #endif
    int T=read();
    while(T--)
    {
        n=read(),m=read();
        tot=ans=0;
        for(int i=1;i<=n;i++) head[i]=siz[i]=vis[i]=0;
        for(int i=1;i<=m;i++) vis[read()]=true;
        for(int i=1;i<n;i++)
        {
            int u=read(),v=read();
            addedge(u,v),addedge(v,u);
        }
        for(int i=1;i<=n;i++)
        {
            if(vis[i])
            {
                dfs(i,0);
                break;
            }
        }
        cout<<ans<<'\n';
    }
    cerr<<endl<<1e3*clock()/CLOCKS_PER_SEC<<"ms";
    return 0;
}