#pragma once
#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <vector>

// For an r-by-b residual panel B, return b-by-b C with B*C=0.
// Its first b-rank(B) columns are a basis of ker(B); later columns are zero.
// Row-major 64-bit words, exactly the worker's coefficient-matrix convention.
namespace residual_cleanup {
struct Result {
    std::vector<std::uint64_t> transform;
    int residual_rank;
    int retained_columns;
};
inline Result kernel_transform(const std::vector<std::uint64_t>& residual, int words) {
    using u64=std::uint64_t;
    if(words<1 || words>16 || residual.size()%words)
        throw std::runtime_error("invalid residual panel dimensions");
    const int b=64*words;
    std::vector<u64> basis(std::size_t(b)*words), value(words);
    std::vector<unsigned char> occupied(b);
    int rank=0;
    for(std::size_t row=0;row<residual.size()/words;++row) {
        std::copy_n(residual.data()+row*words,words,value.data());
        for(;;) {
            int word=0;while(word<words && !value[word])++word;
            if(word==words)break;
            const int pivot=64*word+__builtin_ctzll(value[word]);
            u64* old=basis.data()+std::size_t(pivot)*words;
            if(!occupied[pivot]) {
                std::copy(value.begin(),value.end(),old);occupied[pivot]=1;++rank;break;
            }
            for(int w=word;w<words;++w)value[w]^=old[w];
        }
    }
    Result result{std::vector<u64>(std::size_t(b)*words),rank,b-rank};
    int output=0;
    for(int free=0;free<b;++free)if(!occupied[free]) {
        std::fill(value.begin(),value.end(),0);value[free/64]=u64(1)<<(free%64);
        for(int pivot=b-1;pivot>=0;--pivot)if(occupied[pivot]) {
            const u64* equation=basis.data()+std::size_t(pivot)*words;
            unsigned parity=0;
            for(int w=0;w<words;++w)parity^=__builtin_parityll(equation[w]&value[w]);
            if(parity)value[pivot/64]^=u64(1)<<(pivot%64);
        }
        for(int input=0;input<b;++input)if(value[input/64]>>(input%64)&1)
            result.transform[std::size_t(input)*words+output/64]|=u64(1)<<(output%64);
        ++output;
    }
    if(output!=result.retained_columns)throw std::runtime_error("residual nullity mismatch");
    return result;
}
}
