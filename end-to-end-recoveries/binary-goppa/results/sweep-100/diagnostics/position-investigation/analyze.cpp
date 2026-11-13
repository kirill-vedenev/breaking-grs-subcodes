#define main original_toy_main
#include "../src/toy_goppa.cpp"
#undef main

int main(int argc,char** argv) {
 try {
    require(argc==5,"analyze PUBLIC KERNEL RECOVERED_KEY OUTPUT");
    Public p=read_public(argv[1]);Field f(p.m,p.modulus);
    Binary K=read_kernel(p,argv[2]);Key key=read_key(p,argv[3]);verify_key(p,f,key);
    fs::path output=argv[4];fs::create_directories(output);
    Local local;for(int a=0;a<p.k;++a)for(int b=a+1;b<p.k;++b)local.pairs.emplace_back(a,b);
    int count=1+p.k+int(local.pairs.size()),width=(K->ncols+7)/8;Vec bytes(std::size_t(count)*width);
    require(local_rows(p.k,p.d,p.points[p.held],bytes.data())==0,"Taylor-map construction");
    auto A=binary(count,K->ncols);
    for(int r=0;r<count;++r)for(int c=0;c<K->ncols;++c)
        if((bytes[std::size_t(r)*width+c/8]>>(c%8))&1)mzd_write_bit(A.get(),r,c,1);
    Binary AT(mzd_transpose(nullptr,A.get())),maps(mzd_mul(nullptr,K.get(),AT.get(),0));
    local.J=Matrix(K->nrows,p.k);local.Q=Matrix(K->nrows,int(local.pairs.size()));
    for(int r=0;r<K->nrows;++r) {
        require(!mzd_read_bit(maps.get(),r,0),"nonzero constant Taylor term");
        for(int c=0;c<p.k;++c)local.J(r,c)=mzd_read_bit(maps.get(),r,1+c);
        for(int c=0;c<local.Q.cols;++c)local.Q(r,c)=mzd_read_bit(maps.get(),r,1+p.k+c);
    }
    Matrix Y=p.generator();std::vector<Poly> curve;
    Vec tau(p.k);
    for(int a=0;a<p.k;++a) {
        Vec values(p.n);for(int j=0;j<p.n;++j)values[j]=f.div(Y(a,j),key.lambda[j]);
        Poly poly=interpolate(f,key.support,values);curve.push_back(poly);
        Poly derivative(poly.size());for(int r=1;r<int(poly.size());r+=2)derivative[r-1]=poly[r];
        tau[a]=f.mul(key.lambda[p.held],evaluate(f,derivative,key.support[p.held]));
    }
    Matrix branches(p.m,p.k);for(int s=0;s<p.m;++s)for(int a=0;a<p.k;++a)branches(s,a)=f.power(tau[a],1u<<s);
    Matrix point=from_rows({p.point(p.held)},p.k);
    require(rank(f,stack(point,branches))==p.m+1,"true tangent directions dependent");
    require(multiply(f,local.J,branches.transpose()).zero(),"true tangent gradient mismatch");
    std::vector<std::pair<int,int>> pairs;for(int a=0;a<p.m;++a)for(int b=a;b<p.m;++b)pairs.emplace_back(a,b);
    Matrix substitution(local.Q.cols,int(pairs.size()));
    for(int r=0;r<substitution.rows;++r)for(int c=0;c<substitution.cols;++c) {
        auto [a,b]=local.pairs[r];auto [u,v]=pairs[c];
        substitution(r,c)=f.mul(branches(u,a),branches(v,b));
        if(u!=v)substitution(r,c)^=f.mul(branches(v,a),branches(u,b));
    }
    Matrix restrictions=multiply(f,local.Q,substitution);
    Matrix arcs(p.k,5),monomial_values(K->ncols,5);
    for(int a=0;a<p.k;++a) {
        Poly expansion=compose(f,curve[a],Poly{key.support[p.held],1},5);
        for(int r=0;r<5;++r)arcs(a,r)=f.mul(key.lambda[p.held],coefficient(expansion,r));
    }
    require(series_products(p.k,p.d,p.m,p.modulus,4,arcs.data.data(),monomial_values.data.data())==0,"curve substitution");
    Matrix curve_residual=apply_binary(f,K.get(),monomial_values);
    save_matrix(output/"curve-substitution.txt",curve_residual);
    save_matrix(output/"J.txt",local.J);save_matrix(output/"Q.txt",local.Q);
    save_matrix(output/"branches.txt",branches);save_matrix(output/"branch-quadrics.txt",restrictions);
    save_matrix(output/"branch-quadrics-rref.txt",row_basis(f,restrictions));
    save_matrix(output/"stationary-combinations.txt",nullspace(f,local.J.transpose()));
    Matrix stationary=nullspace(f,local.J.transpose());
    Matrix V(4,4);for(int a=0;a<4;++a)for(int r=0;r<4;++r)V(a,r)=f.power(Byte(a),r);
    Matrix VI=inverse(f,V);
    std::vector<int> cubic_ranks;
    if(rank(f,restrictions)<15)for(int edge=0;edge<6;++edge) {
        int next=(edge+1)%6;Matrix samples(4,K->nrows);
        for(int z=0;z<4;++z) {
            Matrix normalized(p.k,4),products(K->ncols,4);
            for(int a=0;a<p.k;++a) {
                normalized(a,0)=p.point(p.held)[a];
                normalized(a,1)=branches(edge,a)^f.mul(Byte(z),branches(next,a));
                normalized(a,2)=f.sq(normalized(a,1));
            }
            require(series_products(p.k,p.d,p.m,p.modulus,3,normalized.data.data(),products.data.data())==0,"cubic substitution");
            Matrix residual=apply_binary(f,K.get(),products);
            samples.set_row(z,residual.column(3));
        }
        Matrix coefficients=multiply(f,VI,samples).transpose();
        Matrix stationary_coefficients=multiply(f,stationary,coefficients);
        save_matrix(output/("edge-"+std::to_string(edge)+"-cubic.txt"),stationary_coefficients);
        cubic_ranks.push_back(rank(f,stationary_coefficients));
    }
    std::ofstream report(output/"report.json");
    report<<"{\"heldout\":"<<p.held<<",\"support_value\":"<<int(key.support[p.held])
          <<",\"gradient_rank\":"<<rank(f,local.J)<<",\"tangent_span\":"<<rank(f,stack(point,branches))
          <<",\"quadratic_rank\":"<<rank(f,restrictions)<<",\"zero_pairs\":[";
    bool comma=false;
    for(int c=0;c<restrictions.cols;++c)if(pairs[c].first!=pairs[c].second) {
        Vec col=restrictions.column(c);if(std::all_of(col.begin(),col.end(),[](Byte x){return x==0;})){
            if(comma)report<<',';comma=true;report<<'['<<pairs[c].first<<','<<pairs[c].second<<']';
        }
    }
    report<<"],\"continuation_ranks\":[";
    for(int s=0;s<p.m;++s){if(s)report<<',';report<<rank(f,stack(local.J,polar(f,local,branches.row(s))));}
    report<<"],\"curve_nonzero_coefficients_by_order\":[";
    for(int r=0;r<5;++r){if(r)report<<',';Vec col=curve_residual.column(r);report<<std::count_if(col.begin(),col.end(),[](Byte v){return v!=0;});}
    report<<"],\"stationary_cubic_edge_ranks\":[";
    for(std::size_t i=0;i<cubic_ranks.size();++i){if(i)report<<',';report<<cubic_ranks[i];}
    report<<"]}\n";
    std::cout<<"heldout="<<p.held<<" gradient="<<rank(f,local.J)<<" quadratic="<<rank(f,restrictions)<<std::endl;
 } catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}
}
