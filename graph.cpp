#include<bits/stdc++.h>
using namespace std;


/// ////////////////////////////////////////////////////////////////////


using ll = long long;
const int N = 400 , logN = 10 , INF = 1e9;

int n = N;

int low[N] , dfn[N];

int  dfsNumberCounter = 0;
bool UNVISITED        = 0;
int  numSCC           = 0;
int  rootChildren     = 0;

int  dfsRoot;


bool visited[N];
bool articulation_vertex[N];

vector < int > S; // use it like a stack to store SSCs
vector < int > adj[N];// adj list



/// ////////////////////////////////////////////////////////////////////

/// SCC , strongly connected components

void tarjanSCC(int u)
{
    low[u] = dfn[u] = ++ dfsNumberCounter; // low[u] <= dfn[u]
    S.push_back(u); // stores u in a vector based on order of visitation
    visited[u] = 1; // tis vertex is in stack

    for (auto v : adj[u])
    {
        if (dfn[v] == UNVISITED)
            tarjanSCC(v);

        if (visited[v]) // condition for update, v is in stack
            low[u] = min(low[u], low[v]);
    }

    if (low[u] == dfn[u])   // if this is a root (start) of an SCC
    {
        printf("SCC %d:", ++ numSCC); // this part is done after recursion
        while (1)
        {
            int v = S.back();
            S.pop_back(); // pop v out of the stack
            visited[v] = 0; // mark it as outside the stack

            printf(" %d", v);
            if (u == v)
                break;
        }

        printf("\n");
    }
}

/**
    memsetting the values:

    dfn.assign(V, UNVISITED);
    low.assign(V, 0);
    visited.assign(V, 0);
    dfsNumberCounter = numSCC = 0;

    for (int i = 0; i < V; i++)
        if (dfn[i] == UNVISITED)
            tarjanSCC(i);

*/


/// ////////////////////////////////////////////////////////////////////


/// Cut Points & Bridges
void articulationPointAndBridge(int u, int p)
{
    low[u] = dfn[u] = ++ dfsNumberCounter; // low[u] <= dfn[u]

    for (auto v : adj[u])
    {
        if (dfn[v] == UNVISITED)   // a tree edge
        {
            if (u == dfsRoot)
                rootChildren ++; // special case if u is a root

            articulationPointAndBridge(v, u);

            /// cut point
            if (low[v] >= dfn[u]) // for articulation point
                articulation_vertex[u] = true; // store this information first

            /// bridge
            if (low[v] > dfn[u]) // for bridge
                printf(" Edge (%d, %d) is a bridge\n", u, v);

            low[u] = min(low[u], low[v]); // update low[u]
        }
        else if (v != p) // a back edge and not direct cycle
            low[u] = min(low[u], dfn[v]); // update low[u]
    }
}

/**
// inside int main()
    dfsNumberCounter = 0;
    dfn.assign(V, UNVISITED);
    low.assign(V, 0);
    dfs_parent.assign(V, 0);
    articulation_vertex.assign(V, 0);

printf("Bridges:\n");
for (int i = 0; i < V; i++)
    if (dfn[i] == UNVISITED)
    {
        dfsRoot = i;
        rootChildren = 0;
        articulationPointAndBridge(i , 0);/// 0 is a node outside the graph
        articulation_vertex[dfsRoot] = (rootChildren > 1);
    } // special case

printf("Articulation Points:\n");
for (int i = 0; i < V; i++)
    if (articulation_vertex[i])
        printf(" Vertex %d\n", i);

*/

/// ////////////////////////////////////////////////////////////////////


/// LCA

int dep[N];///dep is the depth of a vertex when we run dfs of bfs from the root
int parent[N];


int to [N][logN];
int sum[N][logN];
int mx [N][logN];
int GCD[N][logN];

/// jump k steps above
int jump(int u, int k)
{
    for(int i = 0 ; i < logN ; i ++)
        if(k & (1 << i))
            u = to[u][i];
    return u;
}

/// returns lowest common ancestor of (u , v);
int LCA(int u, int v)
{

    if(u == v)
        return u;

    if(dep[u] < dep[v])
        swap(u, v);

    u = jump(u, dep[u] - dep[v]);
    if(u == v)
        return u;

    for(int i = logN-1 ; i >= 0 ; i --)
        if(to[u][i] != to[v][i])
        {
            u = to[u][i];
            v = to[v][i];
        }

    return to[u][0];

}

/// sum of edges in the path between u ans v through the LCA of u , v
/// same thing to the max_edge or min_edge or GCD
ll SUM(int u , int v){

    ///variable the stores the sum of edges
    ll sumOfEdges = 0ll;
    if(u == v)return sumOfEdges;

    if(dep[u] < dep[v])swap(u , v);

    int k = dep[u] - dep[v];

    // while jumping take the sum with you
    for(int i = 0 ; i < logN ; i ++){
        if(k & (1 << i)){
            sumOfEdges += sum[u][i];
            u           =  to[u][i];
        }
    }

    if(u == v)return sumOfEdges;

    for(int i = logN-1 ; i >= 0 ; i --){
        if(to[u][i] != to[v][i]){
            sumOfEdges += sum[u][i];
            sumOfEdges += sum[v][i];

            u = to[u][i];
            v = to[v][i];
        }
    }

    sumOfEdges += sum[u][0] + sum[v][0];
    return sumOfEdges;
}

/// ////////////////////////////////////////////////////////////////////


/// Sparse Table

/// building sparse table for node , max_edge , sum_of_edges
void build_sparse()
{

    for(int i = 1 ; i <= n ; i ++)
    {
        to[i][0]  = parent[i];
      // sum[i][0]  = //edge_wiegth_between_i_and_parentOf_i;
      //  mx[i][0]  = //edge_wiegth_between_i_and_parentOf_i;
      // GCD[i][0]  = //the gcd between a[i] ans a[ parent[i] ];
    }

    for(int lg = 1 ; lg < logN ; lg ++)
    {
        for(int i = 1 ; i <= n ; i ++)
        {

            to [i][lg] =                         to[ to[i][lg - 1] ][lg - 1];
            sum[i][lg] =    sum[i][lg - 1]    + sum[ to[i][lg - 1] ][lg - 1];
            mx [i][lg] = max(mx[i][lg - 1]    ,  mx[ to[i][lg - 1] ][lg - 1]);
            GCD[i][lg] = __gcd( GCD[i][lg - 1], GCD[ to[i][lg - 1] ][lg - 1] );

        }
    }
}

/// ////////////////////////////////////////////////////////////////////



/// Floyd Warshall’s algorithm applications

int AdjMat[N][N];///for floyd shortest paths
int      p[N][N];///for floyd paths

/// finding all pairs shortest paths
void floyd()
{

// inside int main()
// precondition: AdjMat[i][j] contains the weight of edge (i, j)
// or INF (1B) if there is no such edge
// AdjMat is a 32-bit signed integer array
    for (int k = 1; k <= n; k++) // remember that loop order is k->i->j
        for (int i = 1; i <= n; i++)
            for (int j = 1; j <= n; j++)
                AdjMat[i][j] = min(AdjMat[i][j], AdjMat[i][k] + AdjMat[k][j]);

}

/// are u , v connected ?!
void floyd_connect_or_not()
{

    /// AdjMat is a boolean matrix
    /// if i and j are connected then: AdjMat[i][j] = AdjMat[j][i] = true;
    /// else                           AdjMat[i][j] = AdjMat[j][i] = false;

    for (int k = 1; k <= n; k++)
        for (int i = 1; i <= n; i++)
            for (int j = 1; j <= n; j++)
                AdjMat[i][j] |= (AdjMat[i][k] & AdjMat[k][j]);

}

/// printing the shortest path

void floyd_paths(){

// inside int main()
// let p be a 2D parent matrix, where p[i][j] is the last vertex before j
// on a shortest path from i to j, i.e. i -> ... -> p[i][j] -> j

    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= n; j++)
            p[i][j] = i; // initialize the parent matrix

    for (int k = 1; k <= n; k++)
        for (int i = 1; i <= n; i++)
            for (int j = 1; j <= n; j++) // this time, we need to use if statement
                if (AdjMat[i][k] + AdjMat[k][j] < AdjMat[i][j])
                {
                    AdjMat[i][j] = AdjMat[i][k] + AdjMat[k][j];
                    p[i][j] = p[k][j]; // update the parent matrix
                }
}
//-------------------------------------------------------------------------
// when we need to print the shortest paths, we can call the method below:
void printPath(int i, int j)
{
    if (i != j)
        printPath(i, p[i][j]);
    printf(" %d", j);
}


/// ////////////////////////////////////////////////////////////////////

///dijkstra

int dis[N];

vector < pair < int , int > > adj2[N];

/// ii   = pair < int , int >
/// vii  = vector < ii >
typedef pair < int , int > ii;
typedef vector < ii > vii;
priority_queue < ii , vii , greater < ii > > pq;

void dijkstra(int u){

    // set all distances to INF;
    for(int i = 0 ; i < N ; i ++)dis[i] = INF;

    // set the distance of the starting vertex to 0;
    dis[u] = 0;
    //  {dis , vertex};
    pq.push({0 , u});

    while( !pq.empty() ){

        int u = pq.top().second;
        int d = pq.top().first;
        pq.pop();

        if(d != dis[u])continue; // relax;

        for(auto A : adj2[u]){
            int v = A.second;
            int w = A.first;

            if(w + d < dis[v]){
                dis[v] = d + w;
                pq.push({dis[v] , v});
            }
        }
    }
}

/// ////////////////////////////////////////////////////////////////////


///Topological Sort of a Dag

queue < int > Q;
vector < int > out;
int deg[N] , m;

void  topoSort(){

    // inserting the graph
    for(int i = 0 ; i < m ; i ++){
            int u,v;
            scanf("%d %d",&u,&v);
            adj[u].push_back(v);
            deg[v]++;//increase in_degree of v
    }

    for(int i=1;i<=n;i++)
        if( !deg[i] )//in_degree is 0
            Q.push(i);

        while( !Q.empty() ){
            int u = Q.front();
            Q.pop();

            out.push_back(u);

            for(auto v:adj[u]){
                deg[v]--;//decrease the in_degree of v

                if(!deg[v])//if it became equal to 0 add it to the queue
                    Q.push(v);
            }

        }
    for(int i=0;i<n;i++)
            printf("%d%c",out[i],((i+1==n)?('\n'):(' ')));
}

/// ////////////////////////////////////////////////////////////////////

/// DSU

int pr[N];

int get(int x){
    if(x == pr[x])return x;
    return  pr[x] = get(pr[x]);
}

bool join(int x , int y){
    x = get(x);
    y = get(y);

    if(x == y)return false;

    if(rand() & 1)swap(x , y);

    pr[x] = y;
    return true;
}

/*
in main() {
    for(int i = 0 ; i < N ; i ++)p[i] = i;
}
*/

/// ////////////////////////////////////////////////////////////////////

/// Bellman-Ford
/// Single Source Shortest Paths in graph with negative weight cycles
 
 
vector < pair < int , pair < int , int > > > edges;
 
int Source = 0;
bool angry = false;///detecting a negative cycle in the graph
bool neg[N];/// means this vertex if in negative cycle or not
 
void relax(bool last){
    for(auto E : edges){
        int u = E.second.first;
        int v = E.second.second;
        ll  w = E.first;
 
        if(dis[u] != INF && dis[u] + w < dis[v]){
            dis[v] = dis[u] + w;
            if(last){
                neg[v] = true;
                angry  = true;
            }
        }
    }
}
 
void bellmanFord(){
        dis[Source] = 0;
        for(int i = 1 ; i < n ; i ++){
            relax(0);
        }
 
        relax(1);
}
 
main(){
 
    return 0;
}
