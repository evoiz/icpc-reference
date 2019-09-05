/////////////////////////////////////////////////////
/** segmet tree **/
/**sum ex **/
const int N=1e5+10;
int tree[4*N],a[N];
void build(int p,int l,int r){
    if(l==r){tree[p]=a[l];return;}
    int mid = (l + r) / 2;
    build(p*2,l,mid);
    build(p*2+1,mid+1,r);
    tree[p]=(tree[2*p]+tree[2*p+1]);
}
void update(int p, int l,int r,int i,int x){
    if(l == r){tree[p]=x;return;}
    int mid = (l + r) / 2;
    if( i <= mid)update(p * 2,l , mid, i, x);
    else update(p * 2 + 1, mid + 1, r, i ,x);
    tree[p] =(tree[p * 2]+tree[p * 2 + 1]);
}
int rsq(int p,int l,int r,int i,int j){
    if(l >= i && r <= j){return tree[p];}
    if(r < i || l > j){return 0;}
    int mid = (l + r) / 2;
    int x = rsq(p * 2, l, mid, i, j);
    int y = rsq(p * 2 + 1, mid + 1, r, i, j);
    return x+y;
}
/////////////////////////////////////////////////////
/** segmet tree +lazy **/
/** sum ex**/
const int N=1e5+10;
int tree[4*N],lazy[4*N],n,q,l,r,ty;

void proc(int node, int st, int en) {
if (lazy[node] == 0){return;}
tree[node]+=1ll*(en-st+1)*lazy[node];
if (st != en) {
lazy[node*2]+=1ll*lazy[node];
lazy[node*2+1]+=1ll*lazy[node];
}
lazy[node]=0ll;
}
void update_lazy(int node, int st, int en, int L, int R, ll val) {
proc(node,st,en);
if (st > en || en < L || R < st){return;}
if (L <= st && en <= R) {
lazy[node]+=1ll*val;
proc(node,st,en);
return;
}
update_lazy(node*2, st, (st+en)/2, L, R, val);
update_lazy(node*2+1, (st+en)/2+1, en, L, R, val);
tree[node] =1ll*(tree[node*2]+tree[node*2+1]);
}
ll rsq(int node, int st, int en, int L, int R) {
proc(node, st, en);
if (st > en || en < L || R < st){return 0ll;}
if (L <= st && en <= R){return tree[node];}
ll p1 = rsq(node*2, st, (st+en)/2,L,R);
ll p2 = rsq(node*2+1, (st+en)/2+1,en,L,R);
return 1ll*(p1+p2);
}
/////////////////////////////////////////////////////
struct ProTree {
	vector<int> Tree, Lazy;
	ProTree() {}
	ProTree(int n) {Tree.resize(n*4, 1); Lazy.resize(n*4, 1);}
	
	void Propagate(int node, int st, int en) {
		if (Lazy[node] == 1) return;
		Tree[node] = (1LL * Tree[node] * modPow(Lazy[node], (en - st + 1))) % Mod;
		if (st != en) {
			Lazy[node*2+1] = (1LL * Lazy[node*2+1] * Lazy[node]) % Mod;
			Lazy[node*2+2] = (1LL * Lazy[node*2+2] * Lazy[node]) % Mod;
		}
		Lazy[node] = 1;
	}
	
	void Multiply(int node, int st, int en, int L, int R, int val) {
		Propagate(node, st, en);
		if (st > en || en < L || R < st) return;
		if (L <= st && en <= R) {
			Lazy[node] = (1LL * Lazy[node] * val) % Mod;
			Propagate(node, st, en); return;
		}
		Multiply(node*2+1, st, (st+en)/2+0, L, R, val);
		Multiply(node*2+2, (st+en)/2+1, en, L, R, val);
		Tree[node] = (1LL * Tree[node*2+1] * Tree[node*2+2]) % Mod;
	}
	
	int Product(int node, int st, int en, int L, int R) {
		Propagate(node, st, en);
		if (st > en || en < L || R < st) return 1;
		if (L <= st && en <= R) return Tree[node];
		int p1 = Product(node*2+1, st, (st+en)/2+0, L, R);
		int p2 = Product(node*2+2, (st+en)/2+1, en, L, R);
		return ((1LL * p1 * p2) % Mod);
	}
};

struct OrTree {
	vector<long long> Tree, Lazy;
	OrTree() {}
	OrTree(int n) {Tree.resize(n*4); Lazy.resize(n*4);}
	
	void Propagate(int node, int st, int en) {
		if (Lazy[node] == 0) return;
		Tree[node] |= Lazy[node];
		if (st != en) {
			Lazy[node*2+1] |= Lazy[node];
			Lazy[node*2+2] |= Lazy[node];
		}
		Lazy[node] = 0;
	}
	
	void Update(int node, int st, int en, int L, int R, long long val) {
		Propagate(node, st, en);
		if (st > en || en < L || R < st) return;
		if (L <= st && en <= R) {
			Lazy[node] |= val;
			Propagate(node, st, en); return;
		}
		Update(node*2+1, st, (st+en)/2+0, L, R, val);
		Update(node*2+2, (st+en)/2+1, en, L, R, val);
		Tree[node] = (Tree[node*2+1] | Tree[node*2+2]);
	}
	
	long long Or(int node, int st, int en, int L, int R) {
		Propagate(node, st, en);
		if (st > en || en < L || R < st) return 0;
		if (L <= st && en <= R) return Tree[node];
		long long p1 = Or(node*2+1, st, (st+en)/2+0, L, R);
		long long p2 = Or(node*2+2, (st+en)/2+1, en, L, R);
		return (p1 | p2);
	}

};
/////////////////////////////////////////////////////
/** trie  for char **/
struct Trie{
    Trie *nxt[26];
    Trie *num[14];
    int cnt , isLeaf ;
    Trie(){
        for(int i=0;i<=9;i++){num[i]=0;}
        cnt=isLeaf=0;
        for(int i=0;i<26;i++)
            nxt[i]=0;
        }
};
void add(Trie *tr,int i){
    tr->cnt++;
    if(i==len){tr->isLeaf++;return;}
    if(s[i]>='0' && s[i]<='9'){
       int to = s[i]-'0';
      if(tr->num[to]==0)tr->num[to]=new Trie();
       add(tr->num[to],i+1);
    }
    else{
    int to = s[i]-'a';
    if(tr->nxt[to]==0)tr->nxt[to]=new Trie();
    add(tr->nxt[to],i+1);
    }
}
void del(Trie *tr,int i){
    tr->cnt--;
    if(i==len){tr->isLeaf--;return;}
    if(s[i]>='0' && s[i]<='9'){
      int to = s[i]-'0';
      if(tr->num[to]!=0 && tr->num[to]->cnt>=1){del(tr->num[to],i+1);}
    }
    else{
    int to = s[i]-'a';
    if(tr->nxt[to]!=0 && tr->nxt[to]->cnt>=1){del(tr->nxt[to],i+1);}
    }
}
void get(Trie *tr,int i){
if(i==len){return;}
if(s[i]>='0' && s[i]<='9'){
   int to=s[i]-'0';
   if(tr->num[to]!=0 && tr->num[to]->cnt>=1){
        get(tr->num[to],i+1);
   }
}
else{
int to = s[i]-'a';
if(tr->nxt[to]!=0 && tr->nxt[to]->cnt>=1){
        get(tr->nxt[to],i+1);
}
}
return ;
}
/////////////////////////////////////////////////////
/** trie for bit **/
/** xor ex **/
struct Trie{
    Trie *zero , *one ;
    int cnt;
    Trie(){
        one=zero=NULL;
        cnt=0;
    }
};
void add(Trie *tr,int bt,int x){
    tr->cnt++;
    if(bt<0)return;
    if(x&(1ll<<bt)){
        if( !tr->one )tr->one=new Trie();
        add(tr->one,bt-1,x);
    }
    else{
        if(!tr->zero)tr->zero=new Trie();
        add(tr->zero,bt-1,x);
    }
}
int get(Trie *tr,int i,ll p){
    if( !tr || i<0 ){return 0;}

    if(p&(1ll<<i)){
        if((tr->one)!=0 && (tr->one->cnt>0 ) ) {
            return (1ll<<i)+get(tr->one,i-1,p);
        }
        else{
            return get(tr->zero,i-1,p);
        }
    }
    else{
        if((tr->zero)!=0 && (tr->zero->cnt>0 ) ){
            return get(tr->zero,i-1,p);
        }
        else{
            return (1ll<<i)+get(tr->one,i-1,p);
        }
    }
}
void rem(Trie *tr,int bt,ll x){
    tr->cnt--;
    if(bt<0)return;
    if(x&(1ll<<bt)){
        rem(tr->one,bt-1,x);
    }
    else{
        rem(tr->zero,bt-1,x);
    }
}
/////////////////////////////////////////////////////
/** MO's **/
/** frq of number val in range [l,r] **/ 
const int N=1e6+10;
int fr[N],ans[N],root,L,R,r,l,q,n,t,a[N];

struct node{
int l,r,val,ind;
}qr[N];

bool cmp ( node a , node b ){
if ( a.l/root == b.l/root ) return a.r < b.r ;
return a.l/root < b.l/root ;
}

void add(int ind){fr[a[ind]]++;}
void del(int ind){fr[a[ind]]--;}

int main(){
sc(n);sc(q);
read_arr(a,1,n);
loop(i,1,q){sc(qr[i].l);sc(qr[i].r);sc(qr[i].val);qr[i].ind=i;}
root=ceil(sqrt(n));
sort(qr+1,qr+q+1,cmp);
l=1 ;r=0;
for(int i=1;i<=q;i++){
    L = qr[i].l ,R=qr[i].r ;
    while ( r > R ){del(r);r--;}
    while ( r < R ){r++;add(r);}
    while ( l < L ){del(l);l++;}
    while ( l > L ){l--;add(l);}
    ans[qr[i].ind]=fr[qr[i].val];
}
for(int i=1;i<=q;i++){
 printf("q #%d: %d\n",i,ans[i]);
}
}
/////////////////////////////////////////////////////
/** 2d pr-sum **/
ll query(int x1,int y1,int x2,int y2){
ll ans=sum[x2][y2]-sum[x1-1][y2]-sum[x2][y1-1]+sum[x1-1][y1-1];
return 1ll*ans;
}
void build(){
for(int i=1;i<=n;i++){
   for(int j=1;j<=m;j++)
   sum[i][j]=a[i][j]+sum[i-1][j]+sum[i][j-1]-sum[i-1][j-1];
  }
}
/////////////////////////////////////////////////////

