#define main original_toy_main
#include "../src/toy_goppa.cpp"
#undef main

int main(int argc,char** argv) {
 try {
    require(argc==5,"global_ideal PUBLIC KEY DEGREE OUTPUT");
    Public p=read_public(argv[1]);Field f(p.m,p.modulus);Key key=read_key(p,argv[2]);verify_key(p,f,key);
    int degree=std::stoi(argv[3]);require(degree==3 || degree==4,"degree must be 3 or 4");
    Matrix Y=p.generator();std::vector<Poly> curve;
    for(int a=0;a<p.k;++a) {
        Vec values(p.n);for(int j=0;j<p.n;++j)values[j]=f.div(Y(a,j),key.lambda[j]);
        curve.push_back(interpolate(f,key.support,values));
    }
    std::vector<std::vector<int>> mons;std::vector<int> tmp;
    std::function<void(int)> enumerate=[&](int first) {
        if(int(tmp.size())==degree){mons.push_back(tmp);return;}
        for(int a=first;a<p.k;++a){tmp.push_back(a);enumerate(a);tmp.pop_back();}
    };enumerate(0);int N=mons.size();
    auto equations=binary((degree*p.D()+1)*p.m,N);
    for(int c=0;c<N;++c) {
        Poly value={1};for(int a:mons[c])value=product(f,value,curve[a]);
        for(int r=0;r<int(value.size());++r)for(int b=0;b<p.m;++b)
            if((value[r]>>b)&1)mzd_write_bit(equations.get(),r*p.m+b,c,1);
    }
    int nr=mzd_echelonize_m4ri(equations.get(),1,0);std::vector<int> pivots;std::vector<bool> used(N);
    int cursor=0;for(int r=0;r<nr;++r){while(cursor<N && !mzd_read_bit(equations.get(),r,cursor))++cursor;require(cursor<N,"RREF pivot");pivots.push_back(cursor);used[cursor++]=true;}
    auto K=binary(N-nr,N);int row=0;
    for(int c=0;c<N;++c)if(!used[c]) {
        mzd_write_bit(K.get(),row,c,1);
        for(int r=0;r<nr;++r)if(mzd_read_bit(equations.get(),r,c))mzd_write_bit(K.get(),row,pivots[r],1);
        ++row;
    }
    Binary KT(mzd_transpose(nullptr,K.get())),zero(mzd_mul(nullptr,equations.get(),KT.get(),0));
    require(mzd_is_zero(zero.get()),"global ideal kernel check");
    std::vector<std::pair<int,int>> pairs;for(int a=0;a<6;++a)for(int b=a;b<6;++b)pairs.emplace_back(a,b);
    std::ofstream out(argv[4]);out<<"{\"degree\":"<<degree<<",\"monomials\":"<<N<<",\"binary_constraint_rank\":"<<nr<<",\"ideal_dimension\":"<<N-nr<<",\"positions\":[";
    for(int j=0;j<p.n;++j) {
        Vec point=p.point(j),tau(p.k);Matrix branches(6,p.k);
        for(int a=0;a<p.k;++a) {
            Poly derivative(curve[a].size());for(int r=1;r<int(curve[a].size());r+=2)derivative[r-1]=curve[a][r];
            tau[a]=f.mul(key.lambda[j],evaluate(f,derivative,key.support[j]));
            for(int s=0;s<6;++s)branches(s,a)=f.power(tau[a],1u<<s);
        }
        Matrix values(N,1+p.k+21);
        for(int c=0;c<N;++c) {
            const auto& mon=mons[c];Byte v=1;for(int a:mon)v=f.mul(v,point[a]);values(c,0)=v;
            for(int u=0;u<degree;++u) {
                Byte coefficient=1;for(int w=0;w<degree;++w)if(w!=u)coefficient=f.mul(coefficient,point[mon[w]]);
                values(c,1+mon[u])^=coefficient;
            }
            for(int u=0;u<degree;++u)for(int v=u+1;v<degree;++v) {
                Byte coefficient=1;for(int w=0;w<degree;++w)if(w!=u && w!=v)coefficient=f.mul(coefficient,point[mon[w]]);
                if(!coefficient)continue;
                for(int a=0;a<21;++a) {
                    auto [s,t]=pairs[a];Byte term=f.mul(branches(s,mon[u]),branches(t,mon[v]));
                    if(s!=t)term^=f.mul(branches(t,mon[u]),branches(s,mon[v]));
                    values(c,1+p.k+a)^=f.mul(coefficient,term);
                }
            }
        }
        Matrix image=apply_binary(f,K.get(),values),J(K->nrows,p.k),Q(K->nrows,21);
        for(int r=0;r<K->nrows;++r) {
            require(image(r,0)==0,"global ideal point mismatch");
            for(int a=0;a<p.k;++a)J(r,a)=image(r,1+a);
            for(int a=0;a<21;++a)Q(r,a)=image(r,1+p.k+a);
        }
        if(j)out<<',';out<<"{\"position\":"<<j<<",\"gradient_rank\":"<<rank(f,J)<<",\"quadratic_rank\":"<<rank(f,Q)<<'}';
        std::cout<<"degree="<<degree<<" position="<<j<<" J="<<rank(f,J)<<" Q="<<rank(f,Q)<<std::endl;
    }
    out<<"]}\n";
 } catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}
}
