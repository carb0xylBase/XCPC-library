#include<cstdio>
#include<algorithm>
#include<set>
#include<cmath>
using namespace std;
typedef double db;typedef long long ll;
const int N=2*1e5+10;const db eps=1e-9;
db a[N];db b[N];db r[N];int o[N];int n;db nx;
struct cir//我们这里只存了下标 tp=1/-1 
{
   int tp;int u;cir(int T=0,int N=0){tp=T,u=N;}//这里的计算函数+了一个eps，防止冲突 
   inline db cy(){return b[u]+(db)tp*(sqrt(r[u]*r[u]-(nx-a[u])*(nx-a[u]))+eps);}
   friend bool operator <(cir a,cir b){db jud=a.cy()-b.cy();return jud<-eps;}
};set <cir> s;ll res;
struct pot
{
   int tp;int v;int u;//操作序列点，写了个赋值函数方便压行 
   inline void rv(int T=0,int V=0,int U=0){tp=T;v=V;u=U;}
   friend bool operator <(pot a,pot b){return a.v<b.v;}
}op[2*N];int cnt;
int main()
{freopen("input.txt", "r", stdin);
   scanf("%d",&n);
   for(int i=1;i<=n;i++)//读进来 
   {
   	scanf("%lf%lf%lf",&a[i],&b[i],&r[i]);
   	op[++cnt].rv(1,a[i]-r[i],i);op[++cnt].rv(0,a[i]+r[i],i);
   }sort(op+1,op+cnt+1);//sort一遍 
   for(int i=1;i<=cnt;i++)
   {
   	nx=op[i].v;
   	if(op[i].tp)
   	{
   		int nw=op[i].u;set <cir> ::iterator it;
   		it=s.insert(cir(1,nw)).first;//set的insert返回一个pair类型 
   		if(it==s.begin()){o[nw]=1;}//特判无前驱 
   		else //找到前驱 
   		{
   			--it;if((*it).tp==-1){o[nw]=o[(*it).u]^1;}
   			else {o[nw]=o[(*it).u];}
   		}
   		if(o[nw]){res+=r[nw]*r[nw];}else {res-=r[nw]*r[nw];}//奇加偶减 
   		s.insert(cir(-1,nw));//插入下圆弧 
   	}else{s.erase(cir(1,op[i].u));s.erase(cir(-1,op[i].u));}//大力erase掉即可 
   }printf("%lld",res);return 0;//拜拜程序~ 
}

