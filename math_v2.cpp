#include <bits/stdc++.h>
#define ll long long
using namespace std;
const int mod=1e9+7ll;
//////////////////////////////////////////////////////////////////////
                       /**C(n,k)**/


/** o(n) solution **/
ll comb(int n, int m){
    if(m>n-m){m=n-m;}
    ll ret=1ll;
    //C^{n}_{i} -> C^{n}_{i+1}
    for(int i=0;i<m;++i){
       ret=1ll*ret*(n-i)/(1+i);
   }
    return ret;
}

/** dp solution generate all query in c[N][N] in o(n*m) and query in o(1)**/
/** all pascal tringel**/
const int N=1e2;
ll c[N][N];

void FillLucasTable(){
    memset(c,0,sizeof(c));
    for(int i=0;i<N;++i)c[i][0]=1;
    for(int i=1;i<N;++i)c[i][i]=1;
    for(int i=2;i<N;++i)
        for(int j=1;j<i;++j)
            c[i][j]=(c[i-1][j]+c[i-1][j-1])%mod;
}

/** if query are moded by prime num ex:1e9+7 we can use mod inverse **/
/** mod inverse if we need to find  a/b % mod and gcd(b,mod)==1 then **/
/** a/b % mod= a*pow(b,-1) % mod = a * pow(b, phi(mod) -1 ) % mod **/
/** mod is prime then phi(mod)= mod-1 **/ 
/** c(n,k) % mod = n! /( (n-k)! * (k)! ) % mod = n! * pow( (n-k)!,-1 ) * pow(k!,-1) % mod **/
/** c(n,k) % mod = n! * pow( n-k!,mod-2) * pow(k!,mod-2) % mod **/ 
/** o(n * log(mod) ) **/
const ll mod=1e9+7;
const int N=1e2+10;
ll f[N],inv[N];
ll pow_(ll a,ll b,ll mod=1e9+7){
    ll ans=1ll;
    while(b){
        if(b&1)ans=ans*a%mod;
        a=a*a%mod;
        b>>=1;
    }
    return ans;
}
void init(){
    f[0]=inv[0]=1;
    for(int i=1; i<N; i++){
            f[i]=1ll*i*f[i-1]%mod,inv[i]=pow_(f[i],mod-2);
    }
}
ll C(int n,int m){
    return f[n]*inv[n-m]%mod*inv[m]%mod;
}

//////////////////////////////////////////////////////////////////////
                     /** math law's **/
ll sigma(ll x,ll mod=-1){return 1ll*x*(x+1)/2;}
ll sigma_odd(ll x){return 1ll*((x+1)/2)*((x+1)/2);}
ll sigma_even(ll x){return 1ll*(x/2)*((x/2)+1);}
ll sigma_squares(ll x){return 1ll*x*(x+1)*(2*x+1)/6;}
ll sigma_odd_squares(ll x){x=1ll*(x+1)/2;return 1ll*x*(2*x+1)*(2*x-1)/3;}
ll sigma_even_squares(ll x){x=1ll*x/2;return 1ll*2*x*(x+1)*(2*x+1)/3;}
ll sigma_cubes(ll x){return x*x*(x+1)*(x+1)/4;}
ll sigma_odd_cubes(ll x){x=1ll*(x+1)/2;return 1ll*x*x*(2*x*x-1);}
ll sigma_even_cubes(ll x){x=1ll*x/2;return 1ll*2*x*x*(x+1)*(x+1);}

/**sum of power of number a**/
/** a^0 + a^1 + a^2...+a^b **/
ll sum_of_pow(ll a, ll b){
ll ans=pow_(a,b+1)-1;
ans/=a-1;
return ans;
}
/**with mod && mod inverse**/
ll sum_of_pow(ll a, ll b,ll mod){
ll ans=pow_(a,b+1,mod)-1;
ans=((ans%mod)+mod)%mod;
ans=(ans*pow_(a-1,mod-2,mod))%mod;
return ans;
}

/** negative mod **/
/** ex: (-6) % 4= 2 **/
/** this function can solve positive and negative x but only positive mod **/
ll mod_n(ll x,ll mod){return ((x%mod)+mod)%mod;}

/**first 1 in binary o(1) **/
int msb(unsigned x) {
	union { double a; int b[2]; };
	a = x;
	return (b[1] >> 20) - 1023;
}

/** sum of all digit to n **/
/** ex: sum(11)=(1)+(2)+(3)+(4)+(5)+(6)+(7)+(8)+(9)+(1+0)+(1+1) **/
ll sum_of_digit(ll n){
if(n<10){return n*(n+1)/2;}
int lg=0;
ll tmp=n;
while(tmp/10!=0){lg++;tmp/=10;}
     ll p=ceil(pow(10,lg));
     return ((tmp*45*lg*p/10)+ tmp*(tmp-1)*p/2+ tmp*(n%p+1) + sum_of_digit(n%p));
}

/** gcd between big bumber (saved as string) and number **/
ll big_mod(string s,ll mod){
ll ret=0;
for(int i=0;i<s.length();i++){
  ret*=10;
  ret+=s[i]-'0';
  ret%=mod;
}
return ret;
}

ll big_gcd(string a,ll b){
ll aa=big_mod(a,b);
return __gcd(aa,b);
}

/** check number is fib number **/
/** not very important **/
bool isPerfectSquare(int x){
    int s = sqrt(x);
    return (s*s == x);
}
bool isFibonacci(ll n){
    return isPerfectSquare(5*n*n + 4) ||
           isPerfectSquare(5*n*n - 4);
}

/**NOTE about Fibo:**/
/**The series of last digits repeats with a cycle length of 60**/
/**(Refer this for explanations of this result).**/
/**Every 3-rd Fibonacci number is a multiple of 2**/
/**Every 4-th Fibonacci number is a multiple of 3**/
/**Every 5-th Fibonacci number is a multiple of 5**/
////////////////////////////////////////////////////////////////////////
            /**Euler's totient function**/
                    /** Phi **/
/** phi of number n is: **/
/** the number of k  1<=k<=n where gcd(n,k)=1 **/

/** for example: **/
/** phi(12)=4; **/
/** gcd(12,1)=1  gcd(12,5)=1 gcd(12,7)=1 gcd(12,11)=1 **/
/** so phi(12)={1,5,7,11}=4 **/


/** phi (n)=n*product(1-(1/p) **/
/** where p is all prime number divisor of n **/
/** for example: **/
/** 12=2^2 * 2^3 **/
/** phi(12)=12*(1-(1/2))*(1-(1/3))= 12 *(1/2) * (2/3) **/
/** so we can say phi(n) = n* (p[i]-1) / p[i] : for every prime fact of n p[i] **/

/** O(sqrt (n) ) **/
ll  phi (ll n){
	res=n;
	for (ll i=2;i*i<=n;i++){
		if ( n % i == 0 ){
                while (n%i==0){n/=i;}
               res -= res / i ;
             }
        }
	if(n>1){res-=res/n;}
        return 1ll*res ;
}

/** generate all phi(n) using sieve O(N*log(n)) **/
void sieve_phi() {
   for(int i=1;i<=N;i++){phi[i]=i;}
        for(int i=2;i<N;i++){
           if ( phi[i] == i ) {
              for ( int j=i;j<N;j+=i ) {
                   phi[j]/=1ll*i;
                   phi[j]*=1ll*i-1;
               }
           }
        }
}

/** Application of phi examples **/

/** sum of lcm of n **/
/** sum_of_lcm(n)= sigma( lcm(k,n) ) : 1<=k<=n **/
/** generate all ans in O(N*log N ) **/
/** ans=res[n] **/
/** ex:https://www.codechef.com/problems/SNON02 **/
const int N=1e6+5;
ll phi[N],res[N];
void sum_of_lcm_of_n(){
   sieve_phi();
   for(int i=1;i<N;i++){
        for ( int j=i;j<N;j+=i ){res[j]+=(i*phi[i]);}
    }
   for(int i=1;i<N;i++){res[i]++;res[i]*=i;res[i]/=2ll;}
} 

/** sum of gcd pair of n **/
/** sum_of_gcd_pair(n)= sigma( gcd(i,j) ) : 1<=i,j<=n **/
/** generate all ans in O(N*log N ) **/
/** ans=res[n] **/
/** ex:https://codeforces.com/gym/100145/g**/
void sumOfGcdPairs(){
   sieve_phi();
   for (int i=1;i<N;i++){
       for (int j=2; i*j<N;j++){res[i*j] +=1ll*i*phi[j];}
   }
   for (int i=2; i<N; i++){res[i] += 1ll*res[i-1];}
}
////////////////////////////////////////////////////////////////////////
             /** extended_gcd && solve equations **/

/** ax+by= gcd(a,b)**/
/** function extended_gcd will retrun gcd(a,b) and reference value of &x,&y **/
/** this code can handle negative value **/
ll extended_gcd(ll a,ll b,ll &x,ll &y){
	if(a<0){
        ll ret=extended_gcd(-a,b,x,y);
        x*=-1;
        return ret;
	}
	if(b<0){
        ll ret=extended_gcd(a,-b,x,y);
        y*=-1;
        return ret;
	}
	if(b==0){
		x=1,y=0;
		return a;
	}
	ll g=extended_gcd(b,a%b,y,x);
	y-=(a/b)*x;
	return g;
} 

/** Linear Diophantine equation using extended_gcd**/
/** find solutions for equation ax+by=c**/
/** function give you x0 and y0  you can find other solutions using shift equation **/
/** x=x0+(k*b)/gcd(a,b) **/
/** y=y0+(k*a)/gcd(a,b) **/
bool Linear_Diophantine_equation(ll a,ll b,ll c,ll &x0,ll &y0,ll &g){
	g=extended_gcd(abs(a),abs(b),x0,y0);
	if (c%g!=0){return 0;}
	x0*=c/g ;
	y0*=c/g ;
	if (a<0){x0*=-1ll;}
	if (b<0){y0*=-1ll;}
	return 1;
}

/** shift solution by k **/
/** eqution should by ax+by=gcd(a,b) **/
/** else you shuold div all equation  by gcd(a,b) **/
void shift_solution(ll &x,ll &y,ll a,ll b,ll k){
	x+=k*b ;
	y-=k*a ;
}

/** find number of solutions in range [min_x,max_x] ,[min_y,max_y]**/
ll find_all_solutions(ll a,ll b,ll c,ll minx,ll maxx,ll miny,ll maxy ){
ll x,y,g;
if(!Linear_Diophantine_equation(a,b,c,x,y,g) ){return 0 ;}
a/=g;
b/=g;
ll sign_a=(a>0?1ll:-1ll);
ll sign_b=(b>0?1ll:-1ll);
shift_solution(x,y,a,b,(minx-x)/b);
if(x<minx){shift_solution(x,y,a,b,sign_b);}
if(x>maxx){return 0;}
ll lx1=x;

shift_solution(x,y,a,b,(maxx-x)/b);
if(x>maxx){shift_solution(x,y,a,b,-sign_b);}
ll rx1=x;

shift_solution(x,y,a,b,-(miny-y)/a);
if(y<miny){shift_solution(x,y,a,b,- sign_a);}
if(y>maxy){return 0;}
ll lx2=x;

shift_solution(x,y,a,b,-(maxy-y)/a);
if(y>maxy){shift_solution(x,y,a,b,sign_a);}
ll rx2=x;

if(lx2>rx2){swap(lx2,rx2);}
ll lx=max(lx1,lx2 );
ll rx=min(rx1,rx2 );

return (rx-lx)/abs(b)+1 ;
}

/** for more info **/
/** 
https://translate.google.com/translate?depth=1&hl=en&prev=search&rurl=translate.google.ca&sl=ru&u=http://e-maxx.ru/algo/diofant_2_equation
**/

/** Modular Equation using to solve ax=b (mod m) using extended_gcd**/
/** the function generate all solutions in range 0<=x<mod **/
vector<ll> Modular_Equation(ll a,ll b,ll mod){
set<ll> st;
vector<ll> ret;
ll x,y,g;
g=extended_gcd(a,mod,x,y);
if(b%g!=0){return ret;}
x=((x*b/g)+mod)%mod;
for(int i=0;i<g;i++){
    st.insert((x+i*n/g)%mod);
}
for(auto x:st){ret.pb(x);}
return ret;
}
////////////////////////////////////////////////////////////////////////
                 /**prime nubmer**/

/**seive prime**/
/**generate all prime in range[0,N] **/ 
bool ch[N];
vector<ll>pr;
void sieve(){
for(long long i=2;i<N-3;i++){
    if(!ch[i]){
        pr.pb(i);
        for(long long j=i*2;j<=N-3;j+=i){ch[j]=1;}
    }
 }
}

/**prime factorial for n**/
/**O(sqrt(n)) **/
void f(int n){
    for (int i=2;i*i<=n;i++){
        if (n%i==0){
            cout << i << "^";
            int cnt=0;
            while(n%i==0){
                n/=i;
                cnt++;
            }
            cout << cnt << endl;
        }
    }
    if (n>1){
        cout << n << "^1" << endl;
    }
}

/** when we have a query of prime fact on number x **/
/** we can use sieve to know spf "smallest prime factor" **/
/** then for every query you can div x by spf[x] until x==1 and store prime factor **/
/** generate n*log(n) and query log(x) **/
int spf[N];
void sieve(){
spf[1]=1;
for(int i=2;i<N;i++){spf[i]=i;}
for(int i=4;i<N;i+=2){spf[i]=2;}
for(int i=3;i*i<N;i++){
  if(spf[i]==i){
    for(int j=i*i;j<N;j+=i)
                if(spf[j]==j){spf[j]=i;}
        }
    }
}
vector<int> get(int x) {
vector<int> ret;
while (x!=1){
ret.push_back(spf[x]);
x=x/spf[x];
}
return ret;
}

/** segmentedSieve **/
/** generate prime number between [l,r] in o(n logn) time and (n) memo**/
const int N=1e6+10;
bool ch[N],sg_ch[N];
vector<ll>pr;
void sieve(int n){
pr.clear();
me(ch,0);
for(ll i=2;i<=n;i++){
    if(!ch[i]){
        pr.pb(i);
        for(ll j=i*2;j<=n;j+=i){ch[j]=1;}
    }
 }
}
ll segmentedSieve(ll l,ll h){
	if(l==1){l=2;}
	ll ret=0; 
    ll sq=floor(sqrt(h))+1;  
    sieve(sq); 
    me(sg_ch,0);
     
   for (int i=0;i<(int)pr.size();i++){ 
        ll lo_bl=floor(l/pr[i])*pr[i]; 
        if(lo_bl <  l  ){lo_bl+=pr[i];} 
        if(lo_bl==pr[i]){lo_bl+=pr[i];} 
        for (ll j=lo_bl; j<=h;j+=pr[i]){sg_ch[j-l]=1;} 
    }
    for(ll i=l;i<=h;i++) 
        if(!sg_ch[i-l]){ret++;} 
    return ret;
}
/** miller rabin primality test **/
/** if we need to check prime up to 1e18 **/
/** be careful in 1e18 pow_() function will not wrok because (1e18*1e18) % 1e18 will be overflow **/
/** so we make function called mult_mod(a,b,mod)  retunr a*b % mod with out overflow but complexity will be O( log(n)^2 )**/
ll multMod(ll a,ll b, ll m){
    ll x=0,y=a%m;
    while(b >0){
        if (b%2) x=(x+y)%m;
        y=(y*2)%m;
        b/=2;
    }
    return x;
}
ll pow_(ll a,ll b,ll mod=1e9+7){
    ll ans=1ll;
    a%=mod;
    while(b){
        if(b&1)ans=multMod(ans,a,mod);
        a=multMod(a,a,mod);
        b>>=1;
    }
    return ans;
}
bool miillerTest(ll d, ll n){
    ll a = 2 + rand() % (n - 4);
    ll x = pow_(a, d, n);
    if (x == 1 || x == n-1){return 1;}
    while (d != n-1){
        x =multMod(x,x,n);
        d *= 2;
        if (x == 1){return 0;}
        if (x == n-1){return 1;}
    }
    return 0;
}
bool prime(ll n, ll k){
    if (n <= 1 || n == 4){return 0;}
    if (n <= 3){return 1;}
    ll d = n - 1;
    while (d % 2 == 0){d/= 2;}
    for (ll i = 0; i < k; i++)
         if (!miillerTest(d, n)){return 0;}
    return 1;
}

//////////////////////////////////////////////////////////////////////////////////
                               /**base**/

/** from base 10 to any base **/
/**a0 * base^0 + a1*base^1+....;**/
vector<int> get(int num, int base){
    vector<int> ret;

    if(num==0){ret.push_back(0);return ret;}

    while(num!=0){
        int aux = num%base;
        ret.push_back(aux);
        num /= base;
    }
    return ret;
}

/** div **/
vector<int> div(vector<int> num,int divisor,int base){
reverse(num.begin(),num.end());
vector<int> ret;
int idx=0;
int tmp=num[idx];
while(tmp<divisor){
tmp=tmp*base+(num[++idx]);
}
while (num.size()>idx){
ret.push_back(tmp/divisor);
tmp=(tmp%divisor)*base+num[++idx];
}

if(ret.size()==0){ret.push_back(0);}
reverse(ret.begin(),ret.end());
return ret;
}

/** subtract **/
vector <int> subtract(vector<int > a,vector<int> b,int base){
/**a0 * base^0 + a1*base^1+....;**/
vector<int> ret;
int n=a.size(),m =b.size();
int carry=0;
for(int i=0;i<m;i++){
    int sub = (a[i]-b[i]-carry);
    if(sub<0){sub=sub+base;carry = 1;}
    else{carry=0;}
    ret.push_back(sub);
}
for(int i=m;i<n;i++){
    int sub=(a[i]-carry);
    if(sub<0){sub=sub+base;carry=1;}
    else{carry=0;}
    ret.push_back(sub);
}
return ret;
}

/** add **/
vector<int> add(vector<int > a,vector<int> b,int base){
    /**a0 * base^0 + a1*base^1+....;**/
    vector<int> ret;
    if (a.size() > b.size()){swap(a,b);}
    int n=a.size(),m =b.size();

    int carry = 0;
    for (int i=0;i<n;i++){
        int sum =a[i]+b[i]+carry;
        ret.push_back(sum%base);
        carry =sum/base;
    }
    for(int i=n;i<m;i++){
        int sum = (b[i]+carry);
        ret.push_back(sum%base);
        carry = sum/base;
    }
    if(carry){ret.push_back(carry);}
    return ret;
}

/**left shift**/
vector<int> L_shift(vector<int> a){
vector<int> ret;
ret.pb(0);
for(int i=0;i<a.size();i++){
ret.pb(a[i]);
}
return ret;
}

/** rigth shift **/
vector<int> R_shift(vector<int> a){
a.erase(a.begin());
a.pb(0);
return a;
}
///////////////////////////////////////////////////////////////////////////////////
                     /** the n's Permutation of string S**/
/** given string s of different char what is Permutation number n **/
string solve(int n,string s){
    string ans="";
    int len=s.length();
    vector<char> v;
    for(int i=0;i<len;i++){v.pb(s[i]);}
    ll tmp=len-1,emp=0;
    while(tmp>=1){
    sort(all(v));
    emp=n/fact[tmp];
    n%=fact[tmp];
    ans+=v[emp];
    v.erase(v.begin()+emp);
    tmp-=1ll;
   }
   sort(all(v));
   ans+=v[0];
   return ans;
}
///////////////////////////////////////////////////////////////////////////////////////
                    /** Matrix multiplication algorithm  **/
const int N=2;
ll mod=1e6;
struct matrix{
    int r,c;
    ll mat[N][N];
    matrix(int r_=N,int c_=N){r=r_;c=c_;memset(mat,0,sizeof mat);}
    void print(){
         printf("\n");
            for(int i=0;i<r;i++){
                for(int j=0;j<c;j++){printf("%lld ",mat[i][j]);}
                printf("\n");
        }
        printf("\n");
    }
    matrix operator *(const matrix &B)const{
        matrix ret=matrix(r,B.c);
        for(int i=0;i<r;i++){
            for(int j=0;j<B.c;j++){
                for(int k=0;k<c;k++){
                    ret.mat[i][j]=(1LL*mat[i][k]*B.mat[k][j]+ret.mat[i][j])%mod;
                }
            }
        }
        return ret;
    }
    void init(){
        memset(mat,0,sizeof(mat));
        for(int i=0;i<min(r,c);i++){mat[i][i]=1ll;}
    }
    matrix pow_(ll n){
	   matrix ret;
	   ret.init();
       if(n==0){return ret;}
	   if(n==1){return *this;}
	   matrix P=*this;
	   while(n){
	       if(n&1){ret=ret * P;}
           n=n>>1;
	       P=P*P;
	   }
       return ret;
   }
};
///////////////////////////////////////////////////////////////////////////////////
/* =======> printf <======= */
char ch = 'A';
printf("Character is %c \n", ch);
char str[20] = "fresh2refresh.com";
printf("String is %s \n" , str);
float flt = 10.234;
printf("Float value is %f \n", flt);
double dbl = 20.123456;
printf("Double value is %lf \n", dbl);
int no = 150;
printf("Integer value is %d\n" , no);
printf("Octal value is %o \n", no);
printf("Hexadecimal value is %x \n", no);
printf("%.3lf\n",5+sum/n); //same as setprecision

/* =======> scanf <======= */
char ch;
char str[100];
scanf("%c", &ch);
scanf("%s", &str);
double x,y;
scanf("%lf%lf",&x,&y);
int main(){
    cout<<"Hello ACM!"<< endl;
    cout<<"Math is Math"<<endl;
    cout<<"v.2"<<endl;
    return 0;
}
