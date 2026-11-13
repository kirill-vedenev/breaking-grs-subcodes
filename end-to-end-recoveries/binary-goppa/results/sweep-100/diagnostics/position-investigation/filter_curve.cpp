#define main original_toy_main
#include "../src/toy_goppa.cpp"
#undef main
int main(int argc,char** argv) {
 try {
    require(argc==5,"filter_curve PUBLIC KERNEL KEY OUTPUT");Public p=read_public(argv[1]);Field f(p.m,p.modulus);
    Binary K=read_kernel(p,argv[2]);Key key=read_key(p,argv[3]);verify_key(p,f,key);
    Matrix Y=p.generator(),arcs(p.k,8),values(K->ncols,8);
    for(int a=0;a<p.k;++a) {
        Vec values(p.n);for(int j=0;j<p.n;++j)values[j]=f.div(Y(a,j),key.lambda[j]);
        Poly poly=interpolate(f,key.support,values),local=compose(f,poly,Poly{key.support[p.held],1},8);
        for(int r=0;r<8;++r)arcs(a,r)=f.mul(key.lambda[p.held],coefficient(local,r));
    }
    require(series_products(p.k,p.d,p.m,p.modulus,7,arcs.data.data(),values.data.data())==0,"series products");
    Matrix residual=apply_binary(f,K.get(),values),equations(8*p.m,K->nrows);
    for(int r=0;r<K->nrows;++r)for(int order=0;order<8;++order)for(int b=0;b<p.m;++b)equations(order*p.m+b,r)=(residual(r,order)>>b)&1;
    Matrix combinations=nullspace(f,equations);auto S=binary(combinations.rows,combinations.cols);
    for(int r=0;r<combinations.rows;++r)for(int c=0;c<combinations.cols;++c){require(combinations(r,c)<=1,"nonbinary coefficient");if(combinations(r,c))mzd_write_bit(S.get(),r,c,1);}
    Binary result(mzd_mul(nullptr,S.get(),K.get(),0));std::ofstream out(argv[4],std::ios::binary);out.write("WHK1",4);
    std::uint32_t header[4]={std::uint32_t(p.k),std::uint32_t(p.d),std::uint32_t(K->ncols),std::uint32_t(result->nrows)};
    out.write(reinterpret_cast<char*>(header),sizeof(header));Vec packed((K->ncols+7)/8);
    for(int r=0;r<result->nrows;++r){std::fill(packed.begin(),packed.end(),0);for(int c=0;c<result->ncols;++c)if(mzd_read_bit(result.get(),r,c))packed[c/8]|=1<<(c%8);out.write(reinterpret_cast<char*>(packed.data()),packed.size());}
    require(bool(out),"write failed");
    std::cout<<"{\"input_dimension\":"<<K->nrows<<",\"curve_identity_dimension\":"<<result->nrows<<",\"conditions_added\":"<<rank(f,equations)<<",\"certificate_roots\":256,\"degree_bound\":255}\n";
 } catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}
}
