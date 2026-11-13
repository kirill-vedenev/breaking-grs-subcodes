#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
int main(int argc,char** argv) {
    if(argc!=2)return 2;
    std::ifstream in(argv[1]);int rank,full_n,short_n,held;
    if(!(in>>rank>>full_n>>short_n>>held) || rank>26 || full_n>64 || short_n>64)return 2;
    std::vector<std::uint64_t> full(rank),shortened(rank);
    for(int i=0;i<rank;++i)if(!(in>>full[i]>>shortened[i]))return 2;
    std::vector<std::uint64_t> hf(full_n+1),hs(short_n+1),incident(short_n);
    std::uint64_t a=0,b=0,low=0,held8=0;hf[0]=hs[0]=1;
    for(std::uint64_t v=1;v<(std::uint64_t(1)<<rank);++v) {
        const int bit=__builtin_ctzll(v);a^=full[bit];b^=shortened[bit];
        const int w=__builtin_popcountll(b);++hf[__builtin_popcountll(a)];++hs[w];
        if(w<=12)++low;
        if(w<=8){if(b>>held&1)++held8;for(int j=0;j<short_n;++j)if(b>>j&1)++incident[j];}
    }
    auto list=[](const std::vector<std::uint64_t>& x){std::cout<<'[';for(std::size_t i=0;i<x.size();++i){if(i)std::cout<<',';std::cout<<x[i];}std::cout<<']';};
    std::cout<<"{\"enumerated_dual_vectors\":"<<(std::uint64_t(1)<<rank)<<",\"full_weight_distribution\":";list(hf);
    std::cout<<",\"shortened_weight_distribution\":";list(hs);
    std::cout<<",\"shortened_nonzero_weight_at_most_12\":"<<low<<",\"heldout_incident_weight_at_most_8\":"<<held8;
    std::cout<<",\"coordinate_incidence_weight_at_most_8\":";list(incident);std::cout<<"}\n";
}
