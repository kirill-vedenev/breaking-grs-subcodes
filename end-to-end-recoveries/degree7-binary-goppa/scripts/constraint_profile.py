"""Exact counts and the explicitly heuristic Hilbert-series kernel forecast."""
from math import comb
from algebra import retained_orders, choose


def profile(public):
    k=public['parameters']['k'];d=public['parameters']['degree']
    coefficients=[comb(k,j) for j in range(d+1)]
    active=[];zeros=0;rows=0;nonzeros=0
    for point in public['points']:
        s=point['multiplicity'];weight=point['mask'].bit_count()
        # Unit columns already have multiplicity d-1 for every multilinear h.
        zeros+=d-1 if weight==1 else s
        if weight<d-s+1:continue
        active.append(point)
        exponent=d-s+1
        coefficients=[x-(coefficients[j-exponent] if j>=exponent else 0)
                      for j,x in enumerate(coefficients)]
        levels=retained_orders(d,s)
        rows+=sum(comb(k,u) for u in levels)
        nonzeros+=sum(choose(weight,d-u)*comb(k-d+u,u) for u in levels)
    first=next((j for j,x in enumerate(coefficients) if x<=0),None)
    return dict(active_points=len(active),N=comb(k,d),J=rows,nonzeros=nonzeros,
                forced_zeros=zeros,composition_degree_bound=d*(public['parameters']['n']-2*public['parameters']['t']-1),
                hilbert_coefficients=coefficients,
                first_nonpositive_degree=first,
                predicted_kernel_dimension=coefficients[d] if first is None else 0,
                dimension_forecast_is_heuristic=True,
                dimension_formula='[z^d] ((1+z)^k product_active_points(1-z^(d-s+1)))_+')
