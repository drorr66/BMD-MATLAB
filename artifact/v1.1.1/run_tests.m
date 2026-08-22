function ok = run_tests()
%RUN_TESTS Correctness and canonicality tests for BMDManager v0.8.

fprintf('BMD-MATLAB v0.8 tests\n');
fprintf('======================\n');

nPassed = 0;

%% Constants and monomials
mgr = BMDManager();
assert(isequal(mgr.toDense(mgr.constant(7)), 7));
assert(isequal(mgr.toDense(mgr.monomial(5)), [1 0 0 0 0 0]));
assert(isequal(mgr.toDense(mgr.monomial(6,-3)), [-3 0 0 0 0 0 0]));
pass('constants and monomials');

%% Canonical carry: x*x == x^2
x = mgr.monomial(1);
x2a = mgr.multiply(x,x);
x2b = mgr.monomial(2);
assert(isequal(x2a,x2b), 'Canonical refs differ for x*x and x^2.');
pass('canonical carry x*x == x^2');

%% Paper section 4.2 example: x^5 (x^3 + 2x) = x^8 + 2x^6
F = mgr.monomial(5);
G = mgr.add(mgr.monomial(3), mgr.monomial(1,2));
FG = mgr.multiply(F,G);
assert(isequal(mgr.toDense(FG), [1 0 2 0 0 0 0 0 0]));
pass('paper multiplication example');

%% Different construction paths must give the same canonical reference
mgr2 = BMDManager();
x = mgr2.monomial(1);
one = mgr2.one();
lhs = mgr2.multiply(mgr2.add(x,one), mgr2.subtract(x,one));
rhs = mgr2.subtract(mgr2.monomial(2), one);
assert(isequal(lhs,rhs), 'Canonical refs differ for (x+1)(x-1) and x^2-1.');
pass('canonical equality across construction paths');

%% compact(): remove unreachable construction workspace without changing value
mc = BMDManager();
rc = mc.zero();
for k = 0:100
    rc = mc.add(rc,mc.monomial(k));
end
stBefore = mc.stats(rc);
[mc2,rc2] = mc.compact(rc);
stAfter = mc2.stats(rc2);
assert(isequal(mc.toDense(rc),mc2.toDense(rc2)));
assert(stAfter.total_internal_nodes == stAfter.reachable_nodes);
assert(stAfter.total_internal_nodes == stBefore.reachable_nodes);
assert(stAfter.total_internal_nodes <= stBefore.total_internal_nodes);
pass('compact removes unreachable workspace');

%% v0.4 direct geometric builder: exact, canonical, and garbage-free
for nn = [0 1 2 3 5 10 31 50 100 255]
    md = BMDManager(); rd = md.geometricSum(nn);
    assert(isequal(md.toDense(rd),ones(1,nn+1)));
    sd = md.stats(rd);
    assert(sd.total_internal_nodes == sd.reachable_nodes, ...
        'Direct geometric builder created unreachable nodes for n=%d.',nn);

    mc = BMDManager(); rc = mc.geometricSum(nn);
    rn = mc.zero();
    for kk=0:nn, rn=mc.add(rn,mc.monomial(kk)); end
    assert(mc.same(rc,rn),'Direct and naive geometric builders differ for n=%d.',nn);
end
pass('v0.4 direct geometric builder');

%% v0.4 structural grid: direct sharing equals generic multiplication
for spec = [7 31 8; 8 63 16; 10 127 8]'
    stridePower=spec(1); innerN=spec(2); blocks=spec(3);
    mg=BMDManager();
    direct=mg.geometricGrid(innerN,blocks,stridePower);
    a=mg.geometricSum(innerN);
    b=mg.geometricSumShifted(blocks-1,stridePower+1);
    viaMul=mg.multiply(a,b);
    assert(mg.same(direct,viaMul),'Grid direct/multiply canonical refs differ.');
    stride=2^stridePower;
    expected=zeros(1,(blocks-1)*stride+innerN+1);
    deg=numel(expected)-1;
    exps=[];
    for jj=0:blocks-1, exps=[exps jj*stride+(0:innerN)]; end %#ok<AGROW>
    expected(deg-exps+1)=1;
    assert(isequal(mg.toDense(direct,2e6),expected));
end
pass('v0.4 structural grid sharing');

%% Power-of-two grid has additive (not multiplicative) node growth
for pq = [5 3 7; 6 6 8; 8 8 10]'
    pbits=pq(1); qbits=pq(2); stridePower=pq(3);
    innerN=2^pbits-1; blocks=2^qbits;
    mm=BMDManager(); rr=mm.geometricGrid(innerN,blocks,stridePower); ss=mm.stats(rr);
    assert(ss.reachable_nodes==pbits+qbits, ...
        'Expected p+q structural nodes, got %d for p=%d q=%d.',ss.reachable_nodes,pbits,qbits);
    assert(ss.total_internal_nodes==ss.reachable_nodes);
end
pass('v0.4 additive node growth under block sharing');

%% v0.6 arbitrary support-set builder: canonical and garbage-free
supportSets = {[0 2 5 7],[0 1 4 9 13 31],[3 6 10 15 18 23 30]};
for ii = 1:numel(supportSets)
    ex = supportSets{ii};
    mm = BMDManager();
    direct = mm.indicatorExponents(ex);
    naive = mm.zero();
    for kk = 1:numel(ex)
        naive = mm.add(naive,mm.monomial(ex(kk)));
    end
    assert(mm.same(direct,naive),'indicatorExponents canonical mismatch.');
    sd = mm.stats(direct);
    [mc,rr] = mm.compact(direct);
    sc = mc.stats(rr);
    assert(sc.total_internal_nodes == sc.reachable_nodes);
    assert(isequal(mc.toDense(rr),sparse_terms_to_dense(sparse_terms(uint64(ex),ones(size(ex))))));
end
pass('v0.6 arbitrary exponent-set builder');

%% v0.6 controlled sharing family keeps sparse cardinality fixed
B=8; U=32; K=16; sp=6;
bank = make_sharing_template_bank_v06(B,U,K);
assert(all(sum(bank,2)==K));
assert(size(unique(bank,'rows'),1)==B);
[e1,ids1] = build_sharing_case_v06(bank,1,B,sp);
[e8,ids8] = build_sharing_case_v06(bank,8,B,sp);
assert(numel(e1)==B*K && numel(e8)==B*K);
assert(numel(unique(ids1))==1 && numel(unique(ids8))==8);
mm1=BMDManager(); r1=mm1.indicatorExponents(e1); s1=mm1.stats(r1);
mm8=BMDManager(); r8=mm8.indicatorExponents(e8); s8=mm8.stats(r8);
assert(s8.reachable_nodes > s1.reachable_nodes, ...
    'Removing block-template reuse should increase BMD nodes in the test family.');
assert(isequal(mm1.toDense(r1,1e5),sparse_terms_to_dense(sparse_terms(e1,ones(size(e1))))));
assert(isequal(mm8.toDense(r8,1e5),sparse_terms_to_dense(sparse_terms(e8,ones(size(e8))))));
pass('v0.6 controlled sharing family invariants');


%% v0.7 operation-closure family: fixed input sizes, variable product closure
% P and every Q_s have 256 sparse terms and 8 BMD nodes. Moving Q's bit
% band from full overlap (s=0) to disjoint (s=8) must reduce product/workspace
% complexity while preserving exact multiplication.
closure = struct('shift',{},'result_nodes',{},'new_nodes',{},'sparse_terms',{});
for jj = 1:2
    sp = [0 8]; sp = sp(jj);
    mm = BMDManager();
    pp = mm.geometricSum(255);
    qq = mm.geometricSumShifted(255,sp+1);
    sppStats=mm.stats(pp); sqqStats=mm.stats(qq);
    assert(sppStats.reachable_nodes==8 && sqqStats.reachable_nodes==8);
    [mc,rr] = mm.compact([pp;qq]);
    before = mc.stats();
    prod = mc.multiply(rr(1,:),rr(2,:));
    after = mc.stats(prod);
    spa = build_sparse_family('geometric_sum',255,0);
    spb = build_sparse_shifted_geometric(256,sp);
    spp = sparse_terms_multiply(spa,spb);
    assert(numel(spa.exponents)==256 && numel(spb.exponents)==256);
    for xv=[0.9995 0.9999]
        assert(abs(mc.evaluate(prod,xv)-sparse_terms_evaluate(spp,xv)) <= ...
            1e-8*max(1,abs(sparse_terms_evaluate(spp,xv))));
    end
    closure(jj).shift=sp;
    closure(jj).result_nodes=after.reachable_nodes;
    closure(jj).new_nodes=after.total_internal_nodes-before.total_internal_nodes;
    closure(jj).sparse_terms=numel(spp.exponents);
end
assert(closure(2).result_nodes < closure(1).result_nodes);
assert(closure(2).new_nodes < closure(1).new_nodes);
assert(closure(2).sparse_terms > closure(1).sparse_terms);
pass('v0.7 fixed-input operation-closure family');


%% v0.8 calibrated closure grid: fixed 256x256 inputs, exact new-node targets
pBits=0:7; highSeven=9:15;
overlapBits=[NaN 7 6 5 4 3 2 1 0];
expectedNew=[8 12 16 20 24 28 32 36 40];
expectedResult=[16 16 17 18 19 20 21 22 23];
pEx=build_bitcube_support_v08(pBits);
for ii=1:numel(expectedNew)
    if isnan(overlapBits(ii)), qBits=9:16; else, qBits=[overlapBits(ii) highSeven]; end
    qEx=build_bitcube_support_v08(qBits);
    mm=BMDManager(); pp=mm.indicatorExponents(pEx); qq=mm.indicatorExponents(qEx);
    ps=mm.stats(pp); qs=mm.stats(qq);
    assert(ps.reachable_nodes==8 && qs.reachable_nodes==8);
    [mc,rr]=mm.compact([pp;qq]); before=mc.stats();
    prod=mc.multiply(rr(1,:),rr(2,:)); after=mc.stats(prod);
    assert(before.total_internal_nodes==16);
    assert(after.total_internal_nodes-before.total_internal_nodes==expectedNew(ii));
    assert(after.reachable_nodes==expectedResult(ii));
    if ii==1 || ii==numel(expectedNew)
        spa=sparse_terms(pEx,ones(size(pEx))); spb=sparse_terms(qEx,ones(size(qEx)));
        spp=sparse_terms_multiply(spa,spb);
        assert(isequal(mc.toDense(prod,1e6),sparse_terms_to_dense(spp,1e6)));
    end
end
pass('v0.8 exact closure-threshold calibration grid');

%% compact() with multiple roots preserves canonical identity in one manager
me = BMDManager();
pa = me.zero();
for k = 0:30, pa = me.add(pa,me.monomial(k)); end
q = me.monomial(31);
pb = me.subtract(me.add(pa,q),q);
assert(me.same(pa,pb));
[me2,rr] = me.compact([pa;pb]);
assert(me2.same(rr(1,:),rr(2,:)));
stMulti = me2.stats();
assert(stMulti.total_internal_nodes == me2.reachableNodeCount(rr(1,:)));
pass('multi-root compact preserves canonical equality');

%% Distributivity
mgr3 = BMDManager();
a = mgr3.add(mgr3.monomial(7,2), mgr3.constant(3));
b = mgr3.add(mgr3.monomial(4), mgr3.monomial(1,-2));
c = mgr3.add(mgr3.monomial(3,5), mgr3.constant(-1));
left = mgr3.multiply(a, mgr3.add(b,c));
right = mgr3.add(mgr3.multiply(a,b), mgr3.multiply(a,c));
assert(isequal(left,right), 'Distributivity did not canonicalize to one ref.');
pass('distributivity and canonicalization');

%% Random dense-vs-BMD multiplication checks
rng(1);
for trial = 1:40
    da = randi([0 7]);
    db = randi([0 7]);
    pa = randi([-3 3],1,da+1);
    pb = randi([-3 3],1,db+1);
    if all(pa==0), pa(end)=1; end
    if all(pb==0), pb(end)=1; end
    pa = trimLeading(pa);
    pb = trimLeading(pb);

    m = BMDManager();
    ra = m.fromDense(pa);
    rb = m.fromDense(pb);
    rc = m.multiply(ra,rb);
    got = trimLeading(m.toDense(rc));
    expected = trimLeading(conv(double(pa),double(pb)));
    assert(isequal(got, expected), 'Random multiplication mismatch at trial %d.',trial);
end
pass('40 random multiplications vs conv');

%% Random addition checks
rng(2);
for trial = 1:40
    da = randi([0 10]);
    db = randi([0 10]);
    pa = randi([-5 5],1,da+1);
    pb = randi([-5 5],1,db+1);
    pa = trimLeading(pa); pb = trimLeading(pb);
    n = max(numel(pa),numel(pb));
    expected = [zeros(1,n-numel(pa)) pa] + [zeros(1,n-numel(pb)) pb];
    expected = trimLeading(expected);

    m = BMDManager();
    r = m.add(m.fromDense(pa), m.fromDense(pb));
    got = trimLeading(m.toDense(r));
    assert(isequal(got,double(expected)), 'Random addition mismatch at trial %d.',trial);
end
pass('40 random additions vs dense vectors');


%% Sparse term-list baseline correctness
for nn = [0 1 5 20]
    sg = build_sparse_family('geometric_sum',nn,0);
    assert(isequal(sparse_terms_to_dense(sg),ones(1,nn+1)));
end
sa=sparse_terms(uint64([7 2 0]),[3 -2 5]);
sb=sparse_terms(uint64([5 2 1]),[4 2 -1]);
sadd=sparse_terms_add(sa,sb);
smul=sparse_terms_multiply(sa,sb);
da=sparse_terms_to_dense(sa); db=sparse_terms_to_dense(sb);
N=max(numel(da),numel(db));
expectedAdd=[zeros(1,N-numel(da)) da]+[zeros(1,N-numel(db)) db];
assert(isequal(sparse_terms_to_dense(sadd),trimLeading(expectedAdd)));
assert(isequal(sparse_terms_to_dense(smul),trimLeading(conv(da,db))));
for xv=[-1.25 0 0.5 2]
    assert(abs(sparse_terms_evaluate(sa,xv)-polyval(da,xv))<=1e-12*max(1,abs(polyval(da,xv))));
end
assert(sparse_terms_same(sa,sparse_terms(uint64([7 2 0]),[3 -2 5])));
pass('sparse term-list baseline');

%% Random sparse add/multiply checks vs dense
rng(3);
for trial=1:30
    da=randi([0 8]); db=randi([0 8]);
    pa=randi([-3 3],1,da+1); pb=randi([-3 3],1,db+1);
    pa=trimLeading(pa); pb=trimLeading(pb);
    ea=(numel(pa)-1):-1:0; eb=(numel(pb)-1):-1:0;
    spa=sparse_terms(uint64(ea(pa~=0)),pa(pa~=0));
    spb=sparse_terms(uint64(eb(pb~=0)),pb(pb~=0));
    sadd=sparse_terms_add(spa,spb); smul=sparse_terms_multiply(spa,spb);
    N=max(numel(pa),numel(pb));
    dadd=[zeros(1,N-numel(pa)) pa]+[zeros(1,N-numel(pb)) pb];
    assert(isequal(sparse_terms_to_dense(sadd),trimLeading(dadd)));
    assert(isequal(sparse_terms_to_dense(smul),trimLeading(conv(pa,pb))));
end
pass('30 random sparse operations vs dense');

%% Paper node-count sanity checks (paper counts the terminal as a node)
% Table 1 reports x^n node counts equal to popcount(n)+1.
m = BMDManager();
for n = [50 100 150 200 250 300 350 400 450]
    r = m.monomial(n);
    paperCompatibleCount = m.reachableNodeCount(r) + 1;
    expectedCount = sum(dec2bin(n)=='1') + 1;
    assert(paperCompatibleCount == expectedCount);
end
pass('paper-compatible monomial node counts');

%% More Table 1 structural checks
% sum_{k=0}^n x^k node counts reported in Table 1.
expectedGeo = [11 13 15 15 15 17 17 17 17];
for ii = 1:numel([50 100 150 200 250 300 350 400 450])
    n = [50 100 150 200 250 300 350 400 450];
    nn = n(ii);
    mm = BMDManager(); rr = mm.zero();
    for k = 0:nn
        rr = mm.add(rr,mm.monomial(k));
    end
    assert(mm.reachableNodeCount(rr)+1 == expectedGeo(ii));
end
pass('paper Table 1 geometric-sum node counts');

% Table 1 also reports n+1 nodes for sum k*x^k and (x+1)^n.
for nn = [5 10 20 50]
    mm = BMDManager(); rr = mm.zero();
    for k = 1:nn, rr = mm.add(rr,mm.monomial(k,k)); end
    assert(mm.reachableNodeCount(rr)+1 == nn+1);

    mm = BMDManager(); factor = mm.add(mm.monomial(1),mm.one()); rr = mm.one();
    for k = 1:nn, rr = mm.multiply(rr,factor); end
    assert(mm.reachableNodeCount(rr)+1 == nn+1);
end
pass('paper Table 1 weighted-sum and binomial node counts');

%% Multivariate support and paper Table 2 node-count structure
m = BMDManager();
x1 = m.variable(1);
x2 = m.variable(2);
prod12 = m.multiply(m.add(x1,m.one()), m.add(x2,m.one()));
assert(abs(m.evaluate(prod12,[2 3]) - 12) < 1e-12);
pass('multivariate evaluation');

% Table 2 uses x_k (subscript: distinct variables), not x^k.
% For prod_{k=1}^n (x_k+1)^r the published final *BMD size is r*n+1
% for r=1,4,8. Verify the structural identity on small exact cases.
for rpow = [1 4 8]
    for nv = 1:5
        mm = BMDManager();
        rr = mm.one();
        for vv = 1:nv
            base = mm.add(mm.variable(vv), mm.one());
            rr = mm.multiply(rr, mm.power(base,rpow));
        end
        paperCompatibleCount = mm.reachableNodeCount(rr) + 1;
        assert(paperCompatibleCount == rpow*nv + 1, ...
            'Table 2 node-count mismatch for r=%d,n=%d.',rpow,nv);
    end
end
pass('paper Table 2 multivariate node-count identity');

%% Evaluation
m = BMDManager();
p = [3 -2 0 5 -7];
r = m.fromDense(p);
for xv = [-2 -1 0 0.5 2]
    got = m.evaluate(r,xv);
    expected = polyval(p,xv);
    assert(abs(got-expected) <= 1e-12*max(1,abs(expected)));
end
pass('evaluate vs polyval');

fprintf('\nPASS: %d test groups.\n', nPassed);
ok = true;

    function pass(name)
        nPassed = nPassed + 1;
        fprintf('  PASS  %s\n', name);
    end
end

function p = trimLeading(p)
p = double(p(:).');
idx = find(p~=0,1,'first');
if isempty(idx)
    p = 0;
else
    p = p(idx:end);
end
end
