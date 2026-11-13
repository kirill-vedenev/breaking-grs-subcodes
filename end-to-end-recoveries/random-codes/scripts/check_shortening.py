#!/usr/bin/env python3
"""Public Schur-square checks of a code, its dual, and ordinary shortenings."""
import argparse,hashlib,itertools,json,random,sys,time
from pathlib import Path
from sage.all import GF,PolynomialRing,matrix
ROOT=Path(__file__).resolve().parent

def save(path,data):
 tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(data,default=int,indent=2)+'\n');tmp.replace(path)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('case',type=Path);ap.add_argument('--orderings',type=int,default=8);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 source=args.case/'public.json';data=json.loads(source.read_text());p=data['characteristic'];m=len(data['modulus'])-1
 E=GF(p) if m==1 else GF(p**m,'a',modulus=PolynomialRing(GF(p),'z')(data['modulus']))
 def dec(value):return E(value) if m==1 else sum(E((value//p**i)%p)*E.gen()**i for i in range(m))
 Y=matrix(E,[[dec(a) for a in row] for row in data['matrix']]);H=Y.right_kernel_matrix()
 assert Y.rank()==data['k'] and H.rank()==data['n']-data['k'] and (Y*H.transpose()).is_zero()
 rng=random.Random(9417+data['seed']);started=time.monotonic();records=[]
 for side,code in [('primal',Y),('dual',H)]:
  dimension=code.nrows();n=code.ncols()
  for repetition in range(args.orderings):
   # Reorder randomly, put an information set first, and row-reduce. Every
   # prefix then gives shortening of exactly the displayed codimension.
   permutation=list(range(n));rng.shuffle(permutation)
   permuted=code.matrix_from_columns(permutation)
   pivots=list(permuted.pivots());nonpivots=[j for j in range(n) if j not in pivots]
   ordered=pivots+nonpivots;positions=[permutation[j] for j in ordered]
   G=permuted.matrix_from_columns(ordered).echelon_form()
   assert G.matrix_from_columns(range(dimension)).is_one()
   rows=[]
   for ell in range(dimension):
    A=G[ell:,ell:];r=A.nrows();length=A.ncols();target=min(length,r*(r+1)//2)
    pairs=list(itertools.combinations_with_replacement(range(r),2));rng.shuffle(pairs)
    # Attaining the universal dimension upper bound certifies the full square.
    # Otherwise every product is formed; no sampled plateau is accepted.
    use=min(len(pairs),2*length+16)
    vectors=[A[i].pairwise_product(A[j]) for i,j in pairs[:use]]
    square=matrix(E,vectors);rank=int(square.rank());complete=use==len(pairs)
    if rank<target and not complete:
     square=matrix(E,[A[i].pairwise_product(A[j]) for i,j in pairs]);rank=int(square.rank());complete=True
    rows.append(dict(shortening=ell,length=length,dimension=r,square_rank=rank,random_expected=target,
                     matches_random=rank==target,certificate='attains universal upper bound' if rank==target else 'all products enumerated',
                     products_formed=use if rank==target and not complete else len(pairs)))
   records.append(dict(side=side,ordering=repetition,information_set_original_positions=positions[:dimension],profile=rows))
 result=dict(label=data['label'],seed=data['seed'],parameters={k:data[k] for k in ['q','n','k','D']},public_only=True,
             public_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             orderings_per_side=args.orderings,profiles_checked=sum(len(r['profile']) for r in records),
             all_match_random=all(s['matches_random'] for r in records for s in r['profile']),
             counterexamples=[dict(side=r['side'],ordering=r['ordering'],shortening_positions=r['information_set_original_positions'][:s['shortening']],**s)
                              for r in records for s in r['profile'] if not s['matches_random']],
             sampling_scope='All nonzero dimensions, over the displayed independently shuffled information-set orderings. Not exhaustive over all coordinate subsets.',
             seconds=time.monotonic()-started,records=records)
 args.output.parent.mkdir(parents=True,exist_ok=True);save(args.output,result)
 print(json.dumps({k:v for k,v in result.items() if k!='records'}),flush=True)
if __name__=='__main__':main()
