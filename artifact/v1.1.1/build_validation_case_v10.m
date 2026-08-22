function c = build_validation_case_v10(caseIndex)
%BUILD_VALIDATION_CASE_V10 Deterministic cross-family validation cases.
%
% The v1.0 case set is frozen before MATLAB timing.  It contains regular
% filter/generating-function structures, periodic supports, deterministic
% irregular sparse controls, and weighted coefficient polynomials.

specs = caseSpecs();
if caseIndex < 1 || caseIndex > numel(specs) || caseIndex ~= floor(caseIndex)
    error('BMD:V10Case','caseIndex must be an integer in 1..%d.',numel(specs));
end
s = specs(caseIndex);
c = s;

switch s.kind
    case 'BOX_OVERLAP'
        c.p_sparse = geometricSparse(s.a,1);
        c.q_sparse = geometricSparse(s.b,1);
    case 'BOX_DISJOINT'
        c.p_sparse = geometricSparse(s.a,1);
        c.q_sparse = geometricSparse(s.b,s.step);
    case 'COMB'
        c.p_sparse = geometricSparse(s.a,s.step);
        c.q_sparse = geometricSparse(s.b,s.step2);
    case 'PERIODIC_OVERLAP'
        mask = uint64([0 1 4 5 16 17 20 21]);
        c.p_sparse = repeatMask(mask,s.a,uint64(256));
        c.q_sparse = geometricSparse(s.b,1);
    case 'PERIODIC_DISJOINT'
        mask = uint64([0 1 4 5 16 17 20 21]);
        c.p_sparse = repeatMask(mask,s.a,uint64(256));
        c.q_sparse = geometricSparse(s.b,s.step);
    case 'IRREGULAR'
        ep = deterministicSupport(s.a,s.modulus,s.seed);
        eq = deterministicSupport(s.b,s.modulus,s.seed2);
        c.p_sparse = sparse_terms(ep,ones(1,numel(ep)));
        c.q_sparse = sparse_terms(eq,ones(1,numel(eq)));
    case 'BINOMIAL'
        c.p_sparse = binomialSparse(s.a,1);
        c.q_sparse = binomialSparse(s.b,1);
    case 'DIFFERENCE'
        c.p_sparse = binomialSparse(s.a,-1);
        c.q_sparse = binomialSparse(s.b,-1);
    case 'DENSE_PM1'
        c.p_sparse = densePatternSparse(s.a,s.seed);
        c.q_sparse = densePatternSparse(s.b,s.seed2);
    case 'WEIGHTED_PERIODIC'
        c.p_sparse = weightedPeriodicSparse(s.a,[1 -1 2 -2 3 -3]);
        c.q_sparse = weightedPeriodicSparse(s.b,[2 1 -2 -1]);
    otherwise
        error('BMD:V10Kind','Unknown v1.0 case kind %s.',s.kind);
end

c.p_terms = numel(c.p_sparse.exponents);
c.q_terms = numel(c.q_sparse.exponents);
c.pair_products = double(c.p_terms)*double(c.q_terms);
c.p_degree = degreeOf(c.p_sparse);
c.q_degree = degreeOf(c.q_sparse);
c.result_degree_bound = c.p_degree + c.q_degree;
end

function S=caseSpecs()
S=[ ...
 struct('name','box_overlap_128','family','BOXCAR','application','FIR / moving-average convolution','kind','BOX_OVERLAP','a',128,'b',128,'step',1,'step2',1,'modulus',0,'seed',0,'seed2',0), ...
 struct('name','box_overlap_256','family','BOXCAR','application','FIR / moving-average convolution','kind','BOX_OVERLAP','a',256,'b',256,'step',1,'step2',1,'modulus',0,'seed',0,'seed2',0), ...
 struct('name','box_disjoint_128x320','family','BLOCK_GENERATING','application','blocked generating function / digit-separated kernel','kind','BOX_DISJOINT','a',128,'b',320,'step',256,'step2',1,'modulus',0,'seed',0,'seed2',0), ...
 struct('name','box_disjoint_256x224','family','BLOCK_GENERATING','application','blocked generating function / digit-separated kernel','kind','BOX_DISJOINT','a',256,'b',224,'step',512,'step2',1,'modulus',0,'seed',0,'seed2',0), ...
 struct('name','comb_3x5_128','family','COMB','application','multirate / comb-filter polynomial','kind','COMB','a',128,'b',128,'step',3,'step2',5,'modulus',0,'seed',0,'seed2',0), ...
 struct('name','comb_5x9_192','family','COMB','application','multirate / comb-filter polynomial','kind','COMB','a',192,'b',192,'step',5,'step2',9,'modulus',0,'seed',0,'seed2',0), ...
 struct('name','periodic_mask_overlap','family','PERIODIC_MASK','application','periodic support / repeated spectral mask','kind','PERIODIC_OVERLAP','a',64,'b',128,'step',1,'step2',1,'modulus',0,'seed',0,'seed2',0), ...
 struct('name','periodic_mask_disjoint','family','PERIODIC_MASK','application','periodic support with separated scale','kind','PERIODIC_DISJOINT','a',64,'b',128,'step',32768,'step2',1,'modulus',0,'seed',0,'seed2',0), ...
 struct('name','irregular_sparse_64','family','IRREGULAR_CONTROL','application','irregular sparse control','kind','IRREGULAR','a',64,'b',64,'step',1,'step2',1,'modulus',2048,'seed',17,'seed2',91), ...
 struct('name','irregular_sparse_96','family','IRREGULAR_CONTROL','application','irregular sparse control','kind','IRREGULAR','a',96,'b',96,'step',1,'step2',1,'modulus',4096,'seed',123,'seed2',991), ...
 struct('name','binomial_12x10','family','BINOMIAL','application','binomial transform / repeated first-order factor','kind','BINOMIAL','a',12,'b',10,'step',1,'step2',1,'modulus',0,'seed',0,'seed2',0), ...
 struct('name','binomial_18x14','family','BINOMIAL','application','binomial transform / repeated first-order factor','kind','BINOMIAL','a',18,'b',14,'step',1,'step2',1,'modulus',0,'seed',0,'seed2',0), ...
 struct('name','difference_12x10','family','FINITE_DIFFERENCE','application','finite-difference stencil','kind','DIFFERENCE','a',12,'b',10,'step',1,'step2',1,'modulus',0,'seed',0,'seed2',0), ...
 struct('name','dense_pm1_64','family','DENSE_CONTROL','application','small dense signed-coefficient control','kind','DENSE_PM1','a',64,'b',64,'step',1,'step2',1,'modulus',0,'seed',7,'seed2',19), ...
 struct('name','weighted_periodic_96','family','WEIGHTED_PERIODIC','application','weighted periodic filter kernel','kind','WEIGHTED_PERIODIC','a',96,'b',96,'step',1,'step2',1,'modulus',0,'seed',0,'seed2',0)];
end

function s=geometricSparse(count,step)
if count<1, error('BMD:V10Count','count must be positive.'); end
ex=uint64(0:count-1).*uint64(step);
s=sparse_terms(ex,ones(1,count));
end

function s=repeatMask(mask,blocks,stride)
ex=zeros(1,blocks*numel(mask),'uint64'); pos=0;
for j=0:blocks-1
    idx=pos+(1:numel(mask)); ex(idx)=uint64(j).*stride+mask; pos=pos+numel(mask);
end
s=sparse_terms(ex,ones(1,numel(ex)));
end

function ex=deterministicSupport(count,modulus,seed)
% Small deterministic LCG so the exact support is portable across MATLAB versions.
if count>=modulus, error('BMD:V10Support','count must be smaller than modulus.'); end
ex=zeros(1,count,'uint64'); used=false(1,modulus); x=uint64(seed); k=0;
while k<count
    x=mod(uint64(1664525).*x+uint64(1013904223),uint64(4294967296));
    v=double(mod(x,uint64(modulus)))+1;
    if ~used(v)
        used(v)=true; k=k+1; ex(k)=uint64(v-1);
    end
end
ex=sort(ex,'ascend');
end

function s=binomialSparse(n,signX)
% (1 + signX*x)^n, exact integer coefficients for the small v1.0 cases.
k=0:n; coeff=zeros(1,n+1); coeff(1)=1;
for j=1:n
    coeff(j+1)=coeff(j)*(n-j+1)/j;
end
coeff=coeff.*(signX.^k);
s=sparse_terms(uint64(k),coeff);
end

function s=densePatternSparse(count,seed)
k=0:count-1;
coeff=ones(1,count); coeff(mod(k+seed,3)==0)=-1; coeff(mod(k+2*seed,7)==0)=2;
s=sparse_terms(uint64(k),coeff);
end

function s=weightedPeriodicSparse(count,pattern)
k=0:count-1; coeff=pattern(1+mod(k,numel(pattern)));
s=sparse_terms(uint64(k),coeff);
end

function d=degreeOf(s)
if isempty(s.exponents), d=0; else, d=double(max(s.exponents)); end
end
