// Public-only cubic branch separation and fifth-order lookahead continuation.
// Included after the existing Local structure and polar() helper.
static Vec local_coefficient(const Public& p,const Field& f,const mzd_t* K,
                             const Matrix& arcs,int degree) {
    Matrix products(K->ncols,degree+1);
    require(series_products(p.k,p.d,p.m,p.modulus,degree,arcs.data.data(),products.data.data())==0,"local series products");
    return apply_binary(f,K,products).column(degree);
}
static Vec cubic_obstruction(const Public& p,const Field& f,const mzd_t* K,
                             const Matrix& stationary,const Vec& v) {
    Matrix arcs(p.k,4);Vec point=p.point(p.held);
    for(int a=0;a<p.k;++a){arcs(a,0)=point[a];arcs(a,1)=v[a];arcs(a,2)=f.sq(v[a]);}
    return multiply(f,stationary,local_coefficient(p,f,K,arcs,3));
}
static Matrix fifth_order_matrix(const Public& p,const Field& f,const mzd_t* K,
                                  const Matrix& stationary,const Vec& v) {
    Matrix arcs(p.k,6);Vec point=p.point(p.held);
    for(int a=0;a<p.k;++a){arcs(a,0)=point[a];arcs(a,1)=v[a];arcs(a,2)=f.sq(v[a]);arcs(a,4)=f.power(v[a],4);}
    Vec base=multiply(f,stationary,local_coefficient(p,f,K,arcs,5));
    Matrix result(stationary.rows,p.k);
    for(int a=0;a<p.k;++a) {
        arcs(a,3)=1;Vec shifted=multiply(f,stationary,local_coefficient(p,f,K,arcs,5));arcs(a,3)=0;
        for(int r=0;r<stationary.rows;++r)result(r,a)=base[r]^shifted[r];
    }
    return result;
}
static Local recover_branches_joint(const Public& p,const Field& f,const mzd_t* K) {
    Local local;for(int a=0;a<p.k;++a)for(int b=a+1;b<p.k;++b)local.pairs.emplace_back(a,b);
    int count=1+p.k+int(local.pairs.size()),width=(K->ncols+7)/8;Vec bytes(std::size_t(count)*width);
    require(local_rows(p.k,p.d,p.points[p.held],bytes.data())==0,"Taylor-map construction");
    auto A=binary(count,K->ncols);
    for(int r=0;r<count;++r)for(int c=0;c<K->ncols;++c)
        if((bytes[std::size_t(r)*width+c/8]>>(c%8))&1)mzd_write_bit(A.get(),r,c,1);
    Binary AT(mzd_transpose(nullptr,A.get())),maps(mzd_mul(nullptr,K,AT.get(),0));
    local.J=Matrix(K->nrows,p.k);local.Q=Matrix(K->nrows,int(local.pairs.size()));
    for(int r=0;r<K->nrows;++r) {
        require(!mzd_read_bit(maps.get(),r,0),"holdout polynomial does not vanish at the seed");
        for(int c=0;c<p.k;++c)local.J(r,c)=mzd_read_bit(maps.get(),r,1+c);
        for(int c=0;c<local.Q.cols;++c)local.Q(r,c)=mzd_read_bit(maps.get(),r,1+p.k+c);
    }
    Matrix W=nullspace(f,local.J);require(W.rows==p.m+1,"gradient nullity is not m+1");
    std::vector<Vec> basis={p.point(p.held)};
    for(int r=0;r<W.rows;++r){auto trial=basis;trial.push_back(W.row(r));if(rank(f,from_rows(trial,p.k))>int(basis.size()))basis.push_back(W.row(r));}
    require(int(basis.size())==p.m+1,"gradient quotient basis");basis.erase(basis.begin());local.B=from_rows(basis,p.k);
    local.stationary=nullspace(f,local.J.transpose());
    require(multiply(f,local.stationary,local.J).zero(),"stationary elimination");
    std::vector<std::pair<int,int>> qm;for(int a=0;a<p.m;++a)for(int b=a;b<p.m;++b)qm.emplace_back(a,b);
    Matrix substitution(local.Q.cols,int(qm.size()));
    for(int r=0;r<substitution.rows;++r)for(int c=0;c<substitution.cols;++c) {
        auto [a,b]=local.pairs[r];auto [u,v]=qm[c];substitution(r,c)=f.mul(local.B(u,a),local.B(v,b));
        if(u!=v)substitution(r,c)^=f.mul(local.B(v,a),local.B(u,b));
    }
    Matrix quadrics=multiply(f,local.Q,substitution),quadratic_basis=row_basis(f,quadrics);
    std::vector<std::vector<int>> cm;for(int a=0;a<p.m;++a)for(int b=a;b<p.m;++b)for(int c=b;c<p.m;++c)cm.push_back({a,b,c});
    auto slot=[&](std::vector<int> indices){std::sort(indices.begin(),indices.end());auto it=std::find(cm.begin(),cm.end(),indices);require(it!=cm.end(),"cubic index");return int(it-cm.begin());};
    Matrix cubics(local.stationary.rows,int(cm.size()));
    auto sample=[&](const Vec& coordinates){return cubic_obstruction(p,f,K,local.stationary,multiply(f,local.B.transpose(),coordinates));};
    for(int a=0;a<p.m;++a){Vec e(p.m);e[a]=1;Vec value=sample(e);for(int r=0;r<cubics.rows;++r)cubics(r,slot({a,a,a}))=value[r];}
    Byte beta=2,denominator=f.sq(beta)^beta;
    for(int a=0;a<p.m;++a)for(int b=a+1;b<p.m;++b) {
        Vec e(p.m);e[a]=e[b]=1;Vec u=sample(e);e[b]=beta;Vec w=sample(e);
        for(int r=0;r<cubics.rows;++r) {
            u[r]^=cubics(r,slot({a,a,a}))^cubics(r,slot({b,b,b}));
            w[r]^=cubics(r,slot({a,a,a}))^f.mul(f.power(beta,3),cubics(r,slot({b,b,b})));
            Byte abb=f.div(w[r]^f.mul(beta,u[r]),denominator);
            cubics(r,slot({a,b,b}))=abb;cubics(r,slot({a,a,b}))=u[r]^abb;
        }
    }
    for(int a=0;a<p.m;++a)for(int b=a+1;b<p.m;++b)for(int c=b+1;c<p.m;++c) {
        Vec e(p.m);e[a]=e[b]=e[c]=1;Vec value=sample(e);
        for(int r=0;r<cubics.rows;++r) {
            for(int j:{a,b,c})value[r]^=cubics(r,slot({j,j,j}));
            for(auto pair:std::vector<std::pair<int,int>>{{a,b},{a,c},{b,c}}){auto [u,v]=pair;value[r]^=cubics(r,slot({u,u,v}))^cubics(r,slot({u,v,v}));}
            cubics(r,slot({a,b,c}))=value[r];
        }
    }
    for(Byte c:cubics.data)require(c<=1,"cubic obstruction coefficients must lie in the binary base field");
    Matrix prolonged(quadratic_basis.rows*p.m,int(cm.size()));
    for(int r=0;r<quadratic_basis.rows;++r)for(int j=0;j<p.m;++j)for(int c=0;c<int(qm.size());++c)
        prolonged(r*p.m+j,slot({j,qm[c].first,qm[c].second}))^=quadratic_basis(r,c);
    Matrix ideal=stack(prolonged,cubics),quotient=nullspace(f,ideal);
    require(quotient.rows==p.m,"cubic quotient is not six-dimensional");
    std::vector<Matrix> multiplication;
    for(int a=0;a<p.m;++a) {
        Matrix M(p.m,p.m);for(int j=0;j<p.m;++j)for(int r=0;r<p.m;++r)M(r,j)=quotient(r,slot({0,a,j}));
        multiplication.push_back(M);
    }
    Matrix A1=multiplication[0];require(rank(f,A1)==p.m,"chosen cubic chart misses a branch");
    std::vector<Vec> coordinates;
    for(int mask=2;mask<(1<<p.m) && local.branches.empty();++mask) {
        Matrix A2(p.m,p.m);for(int a=0;a<p.m;++a)if((mask>>a)&1)
            for(std::size_t c=0;c<A2.data.size();++c)A2.data[c]^=multiplication[a].data[c];
        std::vector<Vec> found,found_coordinates;
        for(int theta=0;theta<f.q;++theta) {
            Matrix pencil=A2;for(std::size_t c=0;c<pencil.data.size();++c)pencil.data[c]^=f.mul(Byte(theta),A1.data[c]);
            Matrix left=nullspace(f,pencil.transpose());if(left.rows!=1)continue;
            Vec v=multiply(f,A1.transpose(),left.row(0));found_coordinates.push_back(v);found.push_back(multiply(f,local.B.transpose(),v));
        }
        if(int(found.size())==p.m){local.branches=found;coordinates=found_coordinates;}
    }
    require(int(local.branches.size())==p.m,"could not separate all branches from the joint equations");
    int old_rank=-1;
    for(int s=0;s<p.m;++s) {
        const Vec& v=local.branches[s];const Vec& z=coordinates[s];
        Vec residual=multiply(f,local.J,v);require(std::all_of(residual.begin(),residual.end(),[](Byte x){return !x;}),"branch gradient residual");
        Vec products;for(auto [a,b]:local.pairs)products.push_back(f.mul(v[a],v[b]));residual=multiply(f,local.Q,products);
        require(std::all_of(residual.begin(),residual.end(),[](Byte x){return !x;}),"branch quadratic residual");
        residual=cubic_obstruction(p,f,K,local.stationary,v);require(std::all_of(residual.begin(),residual.end(),[](Byte x){return !x;}),"branch cubic residual");
        Matrix Dq(quadrics.rows,p.m),Dc(cubics.rows,p.m);
        for(int r=0;r<Dq.rows;++r)for(int c=0;c<int(qm.size());++c){auto [a,b]=qm[c];Dq(r,a)^=f.mul(quadrics(r,c),z[b]);Dq(r,b)^=f.mul(quadrics(r,c),z[a]);}
        for(int r=0;r<Dc.rows;++r)for(int c=0;c<int(cm.size());++c)for(int u=0;u<3;++u) {
            Byte value=cubics(r,c);for(int w=0;w<3;++w)if(w!=u)value=f.mul(value,z[cm[c][w]]);Dc(r,cm[c][u])^=value;
        }
        require(rank(f,stack(Dq,Dc))==p.m-1,"joint cubic Jacobian does not restore continuation");
        int r=rank(f,stack(local.J,polar(f,local,v)));if(old_rank<0)old_rank=r;require(r==old_rank,"old ranks differ between Frobenius branches");
    }
    local.extra=fifth_order_matrix(p,f,K,local.stationary,local.branches[0]);
    require(rank(f,stack(stack(local.J,polar(f,local,local.branches[0])),local.extra))==p.k-2,"fifth-order continuation rank");
    local.report.numbers={{"gradient_rank",double(rank(f,local.J))},{"gradient_nullity",double(W.rows)},
        {"quadratic_rank",double(quadratic_basis.rows)},{"cubic_ideal_rank",double(rank(f,ideal))},
        {"cubic_quotient_dimension",double(quotient.rows)},{"branches",double(local.branches.size())},
        {"original_continuation_rank_all_branches",double(old_rank)},{"continuation_rank_all_branches",double(p.k-2)}};
    local.report.checks={{"public_cubic_branch_separation",true},{"stationary_cubics_eliminate_unknown_v3",true},{"fifth_order_lookahead_restores_rank",true}};
    return local;
}

static Jet continue_jet_joint(const Public& p,const Field& f,const mzd_t* K,const Local& local,int depth) {
    require(depth>=2 && depth%2==0,"normalized continuation requests an even final order");
    auto started=Clock::now();Matrix A=stack(stack(local.J,polar(f,local,local.branches[0])),local.extra);FixedSolver solve(f,A);
    require(int(solve.columns.size())==p.k-2,"augmented continuation rank");
    Matrix jets(depth+1,p.k);jets.set_row(0,p.point(p.held));jets.set_row(1,local.branches[0]);
    Vec second=local.branches[0];for(auto& x:second)x=f.sq(x);jets.set_row(2,second);
    std::unique_ptr<void,decltype(&series_cache_free)> cache(series_cache_create(p.k,p.d,p.m,p.modulus,depth+1),series_cache_free);
    require(bool(cache),"series cache allocation");Matrix values(K->ncols,1);
    auto coefficient=[&](int r){require(series_cache_eval(cache.get(),r,values.data.data())==0,"series coefficient evaluation");return apply_binary(f,K,values).column(0);};
    for(int r=0;r<=2;++r) {
        Vec v=jets.row(r);require(series_cache_set(cache.get(),r,v.data())==0,"set initial coefficient");
        auto residual=coefficient(r);require(std::all_of(residual.begin(),residual.end(),[](Byte x){return !x;}),"initial jet substitution");
    }
    int first_full=-1;
    for(int r=3;r<depth;r+=2) {
        Vec even=jets.row((r+1)/2);for(auto& x:even)x=f.sq(x);
        jets.set_row(r+1,even);require(series_cache_set(cache.get(),r+1,even.data())==0,"set known next-even coefficient");
        Vec rhs=coefficient(r),next=coefficient(r+1),extra=multiply(f,local.stationary,coefficient(r+2));
        rhs.insert(rhs.end(),next.begin(),next.end());rhs.insert(rhs.end(),extra.begin(),extra.end());
        Vec odd=solve.solve(rhs);jets.set_row(r,odd);require(series_cache_set(cache.get(),r,odd.data())==0,"update odd coefficient");
        for(int u=r;u<=r+1;++u) {
            auto residual=coefficient(u);require(std::all_of(residual.begin(),residual.end(),[](Byte x){return !x;}),"continued full substitution");
            if(first_full<0 && rank(f,jets.first_rows(u+1))==p.k)first_full=u;
        }
        Vec residual=multiply(f,local.stationary,coefficient(r+2));require(std::all_of(residual.begin(),residual.end(),[](Byte x){return !x;}),"stationary lookahead substitution");
    }
    double recurrence_seconds=elapsed(started);auto tick=Clock::now();Matrix arcs=jets.transpose(),fresh(K->ncols,depth+1);
    require(series_products(p.k,p.d,p.m,p.modulus,depth,arcs.data.data(),fresh.data.data())==0,"fresh full products");
    require(apply_binary(f,K,fresh).zero(),"independent full jet substitution failed");
    for(int r=1;2*r<=depth;++r)for(int a=0;a<p.k;++a)require(jets(2*r,a)==f.sq(jets(r,a)),"binary normalization failed");
    Report report;report.numbers={{"depth",double(depth)},{"continuation_rank",double(solve.columns.size())},
        {"first_full_coefficient_span_order",double(first_full)},{"recurrence_seconds",recurrence_seconds},
        {"fresh_substitution_seconds",elapsed(tick)},{"stationary_lookahead_order",double(depth+1)}};
    report.checks={{"normalized",true},{"full_substitution",true},{"all_polynomials_used",true},{"r_plus_two_stationary_equations",true}};
    return {jets,report};
}
