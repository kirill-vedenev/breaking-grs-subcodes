from pathlib import Path
import json
base=Path(__file__).resolve().parent
old=base.parent/'binary-goppa-e2e/raw'
recovery=old/'diagnostics/seed89-next-heldout-recovery/recovery'
w=(recovery/'reconstructed-curve.txt').read_text().split();r,c=map(int,w[1:3]);a=list(map(int,w[3:]));curve=[a[i*c:(i+1)*c] for i in range(r)]
key=json.loads((recovery/'direct-key.json').read_text())
w=(old/'sweep-100/cases/seed-000089/instance/public.txt').read_text().split();columns=list(map(int,w[9:]));k=28

def mul(a,b):
 c=0
 while b:
  if b&1:c^=a
  b>>=1;a<<=1
  if a&64:a^=67
 return c
M=[[mul(a,b) for b in range(64)] for a in range(64)]
def power(a,n):
 c=1
 while n:
  if n&1:c=M[c][a]
  a=M[a][a];n>>=1
 return c
inv=[0]+[power(a,62) for a in range(1,64)]
def rank(A):
 A=[list(r) for r in A];n=0
 for c in range(len(A[0])):
  p=next((i for i in range(n,len(A)) if A[i][c]),None)
  if p is None:continue
  A[n],A[p]=A[p],A[n];A[n]=[M[inv[A[n][c]]][v] for v in A[n]]
  for i in range(n+1,len(A)):
   if A[i][c]:A[i]=[a^M[A[i][c]][b] for a,b in zip(A[i],A[n])]
  n+=1
  if n==len(A):break
 return n
records=[]
for j,(x,lam) in enumerate(zip(key['support'],key['grs_multiplier'])):
 p=[(columns[j]>>i)&1 for i in range(k)];tau=[]
 for f in curve:
  value=0
  for d in range(len(f)-1,0,-1):value=M[value][x]^(f[d] if d%2 else 0)
  tau.append(M[lam][value])
 radial=next(i for i,v in enumerate(p) if v)
 quotient=[a^M[tau[radial]][b] for a,b in zip(tau,p)]
 pivot=next(i for i,v in enumerate(quotient) if v)
 normalized=[M[inv[quotient[pivot]]][v] for v in quotient]
 fixed={str(d):all(power(v,1<<d)==v for v in normalized) for d in (1,2,3,6)}
 orbit_ranks={str(d):rank([p,tau,[power(v,1<<d) for v in tau]]) for d in (1,2,3)}
 field_degree=next(d for d in (1,2,3,6) if fixed[str(d)])
 assert field_degree==6 and set(orbit_ranks.values())=={3}
 witness=next(v for v in normalized if power(v,4)!=v and power(v,8)!=v)
 records.append({'position':j,'minimal_projective_tangent_field_degree':field_degree,'proper_frobenius_fixed_tests':fixed,'rank_of_point_tangent_and_frobenius_tangent':orbit_ranks,'radial_coordinate':radial,'normalization_coordinate':pivot,'full_degree_coordinate_witness':witness,'witness_fourth_power':power(witness,4),'witness_eighth_power':power(witness,8),'normalized_direction':normalized})
report={'status':'passed','field_modulus':67,'positions':64,'all_projective_tangents_require_GF64':True,'records':records}
(base/'tangent-subfields.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'positions':64,'all_projective_tangents_require_GF64':True,'selected':[records[i] for i in (27,29,42)]},indent=2))
