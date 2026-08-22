function c = build_publication_case_v11(caseIndex)
%BUILD_PUBLICATION_CASE_V11 Frozen publication-validation workload set.
%
% Sixty deterministic polynomial pairs across application-motivated formula
% families plus irregular controls. The set is frozen before MATLAB timing.
S=allSpecs();
if caseIndex<1 || caseIndex>numel(S) || caseIndex~=floor(caseIndex)
    error('BMD:V11Case','caseIndex must be an integer in 1..%d.',numel(S));
end
c=S{caseIndex}; s=c;
switch s.kind
    case 'BOX'
        c.p_sparse=geometricSparse(s.a,1); c.q_sparse=geometricSparse(s.b,1);
    case 'COMB'
        c.p_sparse=geometricSparse(s.a,s.step); c.q_sparse=geometricSparse(s.b,s.step2);
    case 'SEPARATED'
        c.p_sparse=geometricSparse(s.a,1); c.q_sparse=geometricSparse(s.b,s.step);
    case 'PMASK'
        c.p_sparse=repeatMask(s.mask,s.a,uint64(256)); c.q_sparse=geometricSparse(s.b,s.step);
    case 'SUBSET'
        c.p_sparse=subsetProduct(s.weights1); c.q_sparse=subsetProduct(s.weights2);
    case 'BOUNDED'
        c.p_sparse=boundedProduct(s.factors1); c.q_sparse=boundedProduct(s.factors2);
    case 'BINOMIAL'
        c.p_sparse=binomialSparse(s.a,s.sign1); c.q_sparse=binomialSparse(s.b,s.sign2);
    case 'IRREGULAR'
        c.p_sparse=sparse_terms(deterministicSupport(s.a,s.modulus,s.seed),ones(1,s.a));
        c.q_sparse=sparse_terms(deterministicSupport(s.b,s.modulus,s.seed2),ones(1,s.b));
    case 'WEIGHTED'
        c.p_sparse=weightedSparse(s.a,s.pattern1,s.step); c.q_sparse=weightedSparse(s.b,s.pattern2,s.step2);
    otherwise
        error('BMD:V11Kind','Unknown kind %s.',s.kind);
end
c.p_terms=numel(c.p_sparse.exponents); c.q_terms=numel(c.q_sparse.exponents);
c.pair_products=double(c.p_terms)*double(c.q_terms);
c.p_degree=degreeOf(c.p_sparse); c.q_degree=degreeOf(c.q_sparse);
c.result_degree_bound=c.p_degree+c.q_degree;
end

function S=allSpecs()
S={};
% FIR / moving-average convolution: standard boxcar kernels.
for ab=[32 32;64 64;128 64;128 128;256 128;256 256].'
    a=ab(1); b=ab(2); S{end+1}=spec(sprintf('fir_box_%dx%d',a,b),'BOXCAR_FIR','FIR moving-average convolution','APPLICATION_FORMULA','BOX',a,b); %#ok<AGROW>
end
% Comb / CIC-like multirate supports.
C=[32 32 2 3;48 48 2 4;64 64 3 5;64 96 4 8;96 96 5 9;128 64 8 16;128 128 7 11;192 96 16 32];
for k=1:size(C,1)
    a=C(k,1); b=C(k,2); d1=C(k,3); d2=C(k,4); x=spec(sprintf('cic_comb_%dx%d_d%d_d%d',a,b,d1,d2),'CIC_COMB','CIC/comb multirate polynomial','APPLICATION_FORMULA','COMB',a,b); x.step=d1; x.step2=d2; S{end+1}=x; %#ok<AGROW>
end
% Separated-scale/polyphase kernels: low payload times high-rate block index.
C=[64 64 128;64 128 128;128 64 256;128 128 256;128 256 256;256 64 512;256 128 512;256 256 512];
for k=1:size(C,1)
    a=C(k,1); b=C(k,2); st=C(k,3); x=spec(sprintf('polyphase_sep_%dx%d_s%d',a,b,st),'POLYPHASE_SEPARATED','separated-scale polyphase / blocked generating kernel','APPLICATION_FORMULA','SEPARATED',a,b); x.step=st; S{end+1}=x; %#ok<AGROW>
end
% Periodic masks, four overlapping and four separated-scale variants.
M1=uint64([0 1 4 5 16 17 20 21]); M2=uint64([0 2 3 10 11 18 24 25]);
P={M1,16,64,1;M1,32,64,1;M1,32,128,1;M2,32,96,1;M1,16,64,8192;M1,32,64,16384;M1,32,128,32768;M2,32,96,32768};
for k=1:size(P,1)
    mask=P{k,1}; blocks=P{k,2}; b=P{k,3}; st=P{k,4}; x=spec(sprintf('periodic_mask_%02d',k),'PERIODIC_MASK','periodic/repeated spectral mask polynomial','APPLICATION_FORMULA','PMASK',blocks,b); x.mask=mask; x.step=st; S{end+1}=x; %#ok<AGROW>
end
% Subset-sum / 0-1 knapsack generating functions: prod_i (1+x^w_i).
W={ ...
 [1 2 4 8 16 32],[64 128 256 512 1024 2048]; ...
 [1 3 7 15 31 63],[127 255 511 1023 2047]; ...
 [2 5 11 23 47 95],[191 383 767 1535]; ...
 [3 7 13 29 59],[5 11 17 37 73]; ...
 [4 9 19 39 79],[6 13 27 55 111]; ...
 [5 10 20 40 80 160],[7 14 28 56 112 224]; ...
 [1 5 10 25 50],[2 10 20 50 100]; ...
 [6 13 27 54 109],[8 17 35 71 143]};
for k=1:size(W,1)
    x=spec(sprintf('subset_sum_%02d',k),'SUBSET_SUM','0-1 knapsack / subset-sum generating function','APPLICATION_FORMULA','SUBSET',0,0); x.weights1=W{k,1}; x.weights2=W{k,2}; S{end+1}=x; %#ok<AGROW>
end
% Bounded resource/coin-choice generating functions: prod_i sum_{j=0}^m x^(j*d).
F={ ...
 [1 3;4 3;16 3],[64 3;256 3]; ...
 [1 4;5 2;25 2],[125 2;625 2]; ...
 [2 3;7 3;29 2],[5 3;19 2;83 2]; ...
 [3 4;12 3;48 2],[192 2;768 2]; ...
 [1 5;10 3;100 2],[1000 2;10000 1]; ...
 [4 3;20 3;100 2],[6 3;30 3;150 2]};
for k=1:size(F,1)
    x=spec(sprintf('bounded_resource_%02d',k),'BOUNDED_RESOURCE','bounded knapsack/resource generating function','APPLICATION_FORMULA','BOUNDED',0,0); x.factors1=F{k,1}; x.factors2=F{k,2}; S{end+1}=x; %#ok<AGROW>
end
% Binomial / finite-difference transforms.
B=[8 8 1 1;12 10 1 1;16 12 1 -1;18 14 -1 -1;20 8 1 -1;22 10 -1 1];
for k=1:size(B,1)
    a=B(k,1); b=B(k,2); s1=B(k,3); s2=B(k,4); x=spec(sprintf('binomial_diff_%02d',k),'BINOMIAL_DIFFERENCE','binomial / finite-difference transform','APPLICATION_FORMULA','BINOMIAL',a,b); x.sign1=s1; x.sign2=s2; S{end+1}=x; %#ok<AGROW>
end
% Deterministic irregular sparse controls deliberately outside structural families.
I=[24 1024 17 91;32 2048 123 991;40 4096 7 19;48 4096 71 271;56 8192 313 911;64 8192 17 991];
for k=1:size(I,1)
    n=I(k,1); x=spec(sprintf('irregular_control_%02d',k),'IRREGULAR_CONTROL','irregular sparse negative control','CONTROL','IRREGULAR',n,n); x.modulus=I(k,2); x.seed=I(k,3); x.seed2=I(k,4); S{end+1}=x; %#ok<AGROW>
end
% Weighted periodic DSP-style kernels (integer coefficients for exact BMD semantics).
WP={32,[1 -1 2 -2],[2 1 -2 -1],1,1;48,[1 0 -1 2 0 -2],[1 2 -1 -2],1,1;64,[1 -1 1 -1],[2 -2 3 -3],2,3;96,[1 -1 2 -2 3 -3],[2 1 -2 -1],1,1};
for k=1:size(WP,1)
    n=WP{k,1}; x=spec(sprintf('weighted_periodic_%02d',k),'WEIGHTED_PERIODIC','weighted periodic DSP kernel','APPLICATION_FORMULA','WEIGHTED',n,n); x.pattern1=WP{k,2}; x.pattern2=WP{k,3}; x.step=WP{k,4}; x.step2=WP{k,5}; S{end+1}=x; %#ok<AGROW>
end
if numel(S)~=60, error('BMD:V11Specs','Expected exactly 60 frozen cases.'); end
end

function x=spec(name,family,application,sourceBasis,kind,a,b)
x=struct('name',name,'family',family,'application',application,'source_basis',sourceBasis,'kind',kind,'a',a,'b',b, ...
    'step',1,'step2',1,'mask',uint64([]),'weights1',[],'weights2',[],'factors1',[],'factors2',[], ...
    'sign1',1,'sign2',1,'modulus',0,'seed',0,'seed2',0,'pattern1',[],'pattern2',[]);
end
function s=geometricSparse(count,step)
s=sparse_terms(uint64(0:count-1).*uint64(step),ones(1,count));
end
function s=repeatMask(mask,blocks,stride)
ex=zeros(1,blocks*numel(mask),'uint64'); pos=0;
for j=0:blocks-1, idx=pos+(1:numel(mask)); ex(idx)=uint64(j).*stride+mask; pos=pos+numel(mask); end
s=sparse_terms(ex,ones(1,numel(ex)));
end
function s=subsetProduct(weights)
s=sparse_terms(uint64(0),1);
for w=weights, s=localMul(s,sparse_terms(uint64([0 w]),[1 1])); end
end
function s=boundedProduct(factors)
s=sparse_terms(uint64(0),1);
for k=1:size(factors,1)
    d=factors(k,1); m=factors(k,2); f=sparse_terms(uint64(0:m).*uint64(d),ones(1,m+1)); s=localMul(s,f);
end
end
function r=localMul(a,b)
ea=a.exponents; eb=b.exponents; ca=a.coefficients; cb=b.coefficients;
E=ea(:)+eb(:).'; C=ca(:).*cb(:).'; r=sparse_terms(E(:),C(:));
end
function s=binomialSparse(n,signX)
k=0:n; coeff=zeros(1,n+1); coeff(1)=1;
for j=1:n, coeff(j+1)=coeff(j)*(n-j+1)/j; end
coeff=coeff.*(signX.^k); s=sparse_terms(uint64(k),coeff);
end
function ex=deterministicSupport(count,modulus,seed)
if count>=modulus, error('BMD:V11Support','count must be < modulus.'); end
ex=zeros(1,count,'uint64'); used=false(1,modulus); x=uint64(seed); k=0;
while k<count
    x=mod(uint64(1664525).*x+uint64(1013904223),uint64(4294967296)); v=double(mod(x,uint64(modulus)))+1;
    if ~used(v), used(v)=true; k=k+1; ex(k)=uint64(v-1); end
end
ex=sort(ex,'ascend');
end
function s=weightedSparse(count,pat,step)
k=0:count-1; coeff=pat(1+mod(k,numel(pat))); ex=uint64(k).*uint64(step); s=sparse_terms(ex,coeff);
end
function d=degreeOf(s)
if isempty(s.exponents), d=0; else, d=double(max(s.exponents)); end
end
