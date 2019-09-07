#include <bits/stdc++.h>
using namespace std;
typedef long long ll;

const double eps = 1e-10;
const double pi = acos(-1);

struct point{
    double x,y;
    point(){}
    point(double a,double b){ x=a;y=b; }
    void in(){ scanf("%lf %lf",&x,&y); }
    void out(){ printf("{%.9f,%.9f}\n",x,y); }
    double len(){ return hypot(x,y); }
    double len2(){ return x*x+y*y; }
    double angle(){ return atan2(y,x); }
    point operator+(point s){ return point(s.x+x,s.y+y); }
    point operator-(point s){ return point(-s.x+x,-s.y+y); }
    point operator*(double k){ return point(x*k,y*k); }
    point operator/(double k){ return point(x/k,y/k); }
    point  norm(){ return point(x/len(),y/len()); }

};

point vec(point a,point b){ return b-a; }
double dot(point a,point b){ return a.x*b.x + a.y*b.y; }
double crs(point a,point b){ return a.x*b.y - a.y*b.x; }

bool same(point a,point b){
    return fabs( vec(a,b).len() ) <= eps;
}

// compare 2 real numbers a & b
// return +1 if a > b
// return  0 if a = b
// return -1 if a < b
int cmp(double a,double b){
if( fabs(a-b) < eps )return 0;
if( a-b > eps )return 1;
return -1;
}

// compute the distance between c & line defined by (a,b)
double linePointDist(point a,point b,point c){
    return fabs( crs( vec(a,b) , vec(a,c) )/ vec(a,b).len() );
}

double segmentPointDist(point a,point b,point c){
    if( dot(vec(a,b),vec(b,c)) >= 0 )
        return vec(b,c).len();

    if( dot(vec(b,a),vec(a,c)) >= 0 )
        return vec(a,c).len();

    return linePointDist(a,b,c);
}

double polygonArea(vector<point>pol){
double area =0.0;
int N = pol.size();
for(int i=1;i+1<N;i++)
    area += crs( vec(pol[0],pol[i]) , vec(pol[i],pol[i+1]) );
return fabs( area / 2.0 );
}


struct L{
    point p,q;
    double A,B,C;
    L(){}
    L(point aa,point bb){
        p=aa; q=bb;
        A = q.y - p.y;
        B = p.x - q.x;
        C = A * p.x + B * p.y;
    }

    L(double a,double b,double c){
        A=a;
        B=b;
        C=c;
        p=point(1.0,(c-a)/b);
        q=point((c-b)/a,1.0);
        if( same(p,q) )q=point((c-10.0*b)/a,1.0);
    }
// latice points on segment [p,q]
    ll get(){ return 1 + __gcd( (ll)abs( p.x-q.x ) , (ll)abs( p.y-q.y) ); }

// is x in range [l , r]
    bool between(double l,double r,double x){
        if( cmp(l,r)==1 )swap(l,r);
        return cmp(l,x)<=0&&cmp(x,r)<=0;
        }
// is point o on line [p,q]
	bool between(point& o){ return between(p.x,q.x,o.x)&&between(p.y,q.y,o.y); }

    bool intersect(L &o,point &ret){
        double det = A*o.B - o.A*B;
        if( cmp(det,0)==0 )return 0;// parallel or same 2 lines
        double dx = o.B*C - B*o.C;
        double dy = A*o.C - o.A*C;
//        if( dx%det !=0 || dy%det !=0 )return 0;
// not integer coordinates (should be used with int variables det,dx and dy)
        dx /= det;
        dy /= det;
        ret = point(dx,dy);
        return between(ret)&&o.between(ret);
    }
};

void circleFrom3Points(point a,point b,point c,point &cen,double &rad){
    L AB = L(a,b);
    double d1 = - AB.B * ((a.x+b.x)/2.0) + AB.A * ((a.y+b.y)/2.0);
    L t1 = L(-AB.B,AB.A,d1);

    L BC = L(b,c);
    double d2 = - BC.B * ((b.x+c.x)/2.0) + BC.A * ((b.y+c.y)/2.0);
    L t2 = L(-BC.B,BC.A,d2);

    t1.intersect(t2,cen);
    rad = vec(cen,a).len();
}

point rot(point p,double t){
    return point( p.x*cos(t) - p.y*sin(t) , p.x*sin(t) + p.y*cos(t) );
}

point rotAboutPoint(point p,double t,point q){
    return q+rot(p-q,t);
}

bool cmp1(point a,point b){
    return cmp(a.x,b.x)==-1 ||( cmp(a.x,b.x)==0 && cmp(a.y,b.y)==-1 );
}

vector<point> convexHull(vector<point>pol){
sort(pol.begin(),pol.end(),cmp1);
// remove duplicate points
vector<point>tmp;
tmp.push_back(pol[0]);
int N = pol.size();
for(int i=1;i<N;i++)
    if( !same(pol[i-1],pol[i]) ){tmp.push_back(pol[i]);}
pol=tmp;

N=pol.size();// again
vector<point>up,dn;
for(int i=0;i<N;i++){
        while( dn.size() >=2 && crs(vec(dn[dn.size()-2],dn[dn.size()-1]),vec(dn[dn.size()-1],pol[i])) <0 ){dn.pop_back();}
        dn.push_back(pol[i]);
}

for(int i=N-1;i>=0;i--){
        while( up.size() >=2 && crs(vec(up[up.size()-2],up[up.size()-1]),vec(up[up.size()-1],pol[i])) <0 )
            up.pop_back();
    up.push_back(pol[i]);
    }

dn.pop_back();
up.pop_back();

vector<point>cvx;
for(auto p:dn)cvx.push_back(p);
for(auto p:up)cvx.push_back(p);
return cvx;
}

bool pointInPolygon(vector<point>&pol,point &a){
int N = pol.size();
double tot =0.0;
for(int i=0;i<N;i++)
    tot += vec(a,pol[(i+1)%N]).angle() - vec(a,pol[i]).angle();

return (cmp(tot,eps) != 0);
}

bool pointInTriangle(point a,point b,point c,point p){
    return crs( vec(a,b) , vec(a,p) ) >=0
        && crs( vec(b,c) , vec(b,p) ) >=0
        && crs( vec(c,a) , vec(c,p) ) >=0 ;
}

bool pointInConvex(vector<point>&pol,point &a){
int N = pol.size();
if( crs( vec(pol[0],pol[1]) , vec(pol[0],a) ) < 0 )return 0;
if( crs( vec(pol[0],pol[N-1]) , vec(pol[0],a) ) > 0 )return 0;
int low =1,high = N-2;
while(low < high){
    int mid = (low+high+1)>>1;
    if( crs(vec(pol[0],pol[mid]),vec(pol[0],a)) >0 )low=mid;
    else high=mid-1;
}
return pointInTriangle(pol[0],pol[low],pol[low+1],a);
}
//////////////////////////////////////////////////////////////////////////////
/**tringel**/
double AreaHeron(double const &a, double const &b, double const &c){
    double s=(a+b+c)/2;
    return sqrt(s*(s-a)*(s-b)*(s-c));
}

double Circumradius(const double &a, const double &b, const double &c){
    return a*b*c/4/AreaHeron(a,b,c);
}

double Circumradius(const point &P1, const point &P2, const point &P3){
    return (P2-P1).abs()*(P3-P1).abs()*(P3-P2).abs()/4/fabs(signed_area(P1,P2,P3));
}

double Inradius(const double &a, const double &b, const double &c){
    return 2*AreaHeron(a,b,c)/(a+b+c);
}
//////////////////////////////////////////////////////////////////////////////////
/**points**/
const double eps = 1e-9;

struct point{
    double x,y;
    
    point(){}
    point(double _x, double _y) : x(_x), y(_y){}
    
    double cross(point P){
        return x * P.y - y * P.x;
    }

    bool operator < (const point &p) const{
        if(fabs(x-p.x)>eps) return x<p.x;
        return y>p.y;
    }
};

double cross(point a, point b){
    return a.x * b.y - a.y * b.x;
}

bool polar_cmp(point a, point b){
    if(a.x >= 0 && b.x < 0) return true;
    if(a.x < 0 && b.x >= 0) return false;
    if(a.x == 0 && b.x == 0){
        if(a.y > 0 && b.y < 0) return false;
        if(a.y < 0 && b.y > 0) return true;
    }
    return cross(a,b) > 0;
}

struct line{
    point p1,p2;
    
    line(){}
    
    line(point _p1, point _p2){
        p1 = _p1; p2 = _p2;
        if(p1.x>p2.x) swap(p1,p2);
    }
};

double PointToLineDist(const point &P, const line &L){
    return 2 * fabs(signed_area(L.p1,L.p2,P)) / (L.p2 - L.p1).abs();
}
double PointToSegmentDist(const point &P, const line &L){
    point v = L.p2 - L.p1,w = P - L.p1;
    
    double aux1 = w.dot(v);
    if(aux1 < eps) return (P-L.p1).abs();
    
    double aux2 = v.dot(v);
    if(aux2 <= aux1+eps) return (P-L.p2).abs();
    
    return PointToLineDist(P,L);
}

bool cross(point a, point b, point c, point d){
    return area(a,b,c) * area(a,b,d) < eps && area(c,d,a) * area(c,d,b) < eps;
}

point intersection(point a, point b, point c, point d){
    return a + (b - a) * ((c - a).cross(d - c) / (b - a).cross(d - c));
}

point projection(point x, point a, point b){
    return a + (b - a) * ((x - a).dot(b - a) / (b - a).dot(b - a));
}

point reflection(point x, point a, point b){
    return x + (projection(x,a,b) - x) * 2;
}

//verdadero : sentido anti-horario, Complejidad : O(n)
bool ccw(const vector<point> &poly){
    //primero hallamos el punto inferior ubicado mas a la derecha
    int ind = 0,n = poly.size();
    double x = poly[0].x,y = poly[0].y;

    for(int i = 1;i < n;i++){
        if(poly[i].y > y) continue;
        if(fabs(poly[i].y - y) < eps && poly[i].x< x) continue;
        ind = i;
        x = poly[i].x;
        y = poly[i].y;
    }
    
    if(ind == 0) return ccw(poly[n - 1],poly[0],poly[1]);
    return ccw(poly[ind-1],poly[ind],poly[(ind + 1) % n]);
}

bool isInConvex(vector <point> &A, const point &P){
    // polygon should be in ccw order
    // if(!ccw(A)) reverse(A.begin(),A.end());
    int n = A.size(),lo = 1,hi = A.size() - 1;
    
    if(area(A[0],A[1],P) <= 0) return 0;
    if(area(A[n - 1],A[0],P) <= 0) return 0;
    
    while(hi - lo > 1){
        int mid = (lo + hi) / 2;
        
        if(area(A[0],A[mid],P) > 0) lo = mid;
        else hi = mid;
    }
    
    return area(A[lo],A[hi],P) > 0;
}



bool PointInsidePolygon(const point &P, const vector<point> &poly){
    int n = poly.size();
    bool in = 0;
    
    for(int i = 0,j = n - 1;i < n;j = i++){
        double dx = poly[j].x - poly[i].x;
        double dy = poly[j].y - poly[i].y;
        
        if((poly[i].y <= P.y + eps && P.y < poly[j].y) ||
            (poly[j].y <= P.y + eps && P.y < poly[i].y))
            if(P.x - eps < dx * (P.y-poly[i].y) / dy+poly[i].x)
                in ^= 1;
    }
    
    return in;
}

//valor positivo : vertices orientados en sentido antihorario
//valor negativo : vertices orientados en sentido horario
double signed_area(const vector<point> &poly){
  int n = poly.size();
  if(n < 3) return 0.0;
  
  double S = 0.0;
  
  for(int i = 1;i <= n;++i)
    S += poly[i % n].x * (poly[(i + 1) % n].y - poly[i - 1].y);
  
  return S / 2;
}

