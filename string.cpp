/** KMP **/
/** find all pos of pattern P in sting s **/
/** O(p.length() + s.length() )**/
vector<int> prefixFunction(string p){
    int n = p.size(), k = 0;
    vector<int> ret(n);
    ret[0]=0;
    for(int i=1;i<n;++i){
        while(k>0 && p[k]!=p[i]){k=ret[k-1];}
        if(p[k]==p[i]){++k;}
        ret[i]=k;
    }
    return ret;
}
vector<int> KMP(string p, string s){
    int m=p.size(),n= s.size(),k=0;
    vector<int> ret;
    vector<int> pre=prefixFunction(p);
    for(int i=0;i<n;++i){
        while(k>0 && p[k]!=s[i]){k=pre[k-1];}
        if(p[k]==s[i]){++k;}
        if(k==m){
            ret.push_back(i-k+1);
            k=pre[k-1];
        }
    }
    return ret;
}
////////////////////////////////////////////////////////////////////////
