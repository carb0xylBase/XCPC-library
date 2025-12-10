#include<bits/stdc++.h>
using namespace std;
using ll=long long;
using ull=unsigned long long;
const int N=400020;
const int inf=1e9;

int mod,a[N],sum[N],sum2[N],sum3[N];
inline int ksm(int x,int k) {int s=1;while(k){if(k&1) s=1ll*s*x%mod;x=1ll*x*x%mod;k>>=1;}return s;}
void sol()
{
    int n,m;
    cin>>n>>m>>mod; 
    int invm=ksm(m,mod-2),now=1;
    int ans=0;
    for(int i=1;i<=n;i++)
    {
        now=1ll*now*invm%mod;
        sum[i]=(sum[i-1]+now)%mod;  
        sum2[i]=(sum2[i-1]+1ll*now*(i-1)%mod)%mod;
        sum3[i]=(sum3[i-1]+1ll*now*i%mod)%mod;
    }
    int len=2*n-1,tot=0;
    for(int i=1;i<=len;i++)
    {
        int k=min(i/2,(len-i+1)/2);
        a[i]=(sum[k]+(i%2==1))%mod;
        tot=(tot+a[i])%mod;

        if(i%2==1) ans=(ans+2ll*sum3[k]%mod)%mod;
        else ans=(ans+2ll*sum2[k]%mod)%mod;
        ans=(ans+sum[k]+(i%2==1))%mod;
    }
    for(int i=1;i<=len;i++)
        ans=(ans+1ll*(tot-a[i]+mod)%mod*a[i]%mod)%mod;
    cout<<ans<<endl;
}
signed main()
{
    freopen("input.txt", "r", stdin);
    ios::sync_with_stdio(0);
    cin.tie(0);cout.tie(0);

    int T;cin>>T;
    while(T--)
        sol();
    return 0;
}