%% Polyphase FIR application validation for exact DBB closure
% External to the frozen v1.1.1 60-case benchmark.
% Source design: MathWorks FIR Decimation for FPGA example, decimation by 8.

clear; clc;
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
frozenDir = fullfile(repoRoot,'artifact','v1.1.1');
addpath(frozenDir);

outDir = fullfile(fileparts(mfilename('fullpath')),'results');
if ~exist(outDir,'dir'), mkdir(outDir); end

%% Fixed application configuration (do not tune after seeing results)
M = 8;
s = 3;
qBits = 15;
coeffs = firpm(30,[0 0.1 0.2 0.5]*2,[1 1 0 0]);
q = round(coeffs * 2^qBits);
if all(q==0), error('Application:FIRQuantization','Q15 quantization produced all zeros.'); end
P = sparse_terms(0:numel(q)-1,q);

%% A. Search the complete quantized FIR for nontrivial exact DBB factorization
searchRows = [];
for ss = 1:4
    stride = 2^ss;
    [ok,Fsp,Gsp,why] = exact_dbb_factorization(P,stride);
    if ok
        [mgr,refs] = build_bmd_pair_from_sparse_v10(Fsp,Gsp);
        nF = mgr.reachableNodeCount(refs(1,:));
        nG = mgr.reachableNodeCount(refs(2,:));
        st0 = mgr.stats();
        prod = mgr.multiply(refs(1,:),refs(2,:));
        st1 = mgr.stats(prod);
        nNew = st1.nodes_created - st0.nodes_created;
        nProd = mgr.reachableNodeCount(prod);
        sparseProd = sparse_terms_multiply(Fsp,Gsp);
        exactOK = sparse_terms_same(sparseProd,P);
        predNewOK = (nNew == nF);
        predSizeOK = (nProd == nF+nG);
        nontrivial = numel(Fsp.exponents)>=2 && nF>=2;
    else
        nF=NaN; nG=NaN; nNew=NaN; nProd=NaN;
        exactOK=false; predNewOK=false; predSizeOK=false; nontrivial=false;
    end
    searchRows = [searchRows; table(ss,stride,ok,nontrivial,nF,nG,nNew,nProd, ...
        predNewOK,predSizeOK,exactOK,string(why), ...
        'VariableNames',{'s','stride','dbb_factorization','nontrivial_low_factor', ...
        'V_F','V_G','N_new','V_product','prediction_new_ok','prediction_size_ok', ...
        'exact_product_ok','note'})]; %#ok<AGROW>
end
writetable(searchRows,fullfile(outDir,'whole_filter_dbb_search.csv'));

%% B. Natural M=8 polyphase branches x^k H_k(x^8)
branchRows = [];
for k = 0:M-1
    idx = (k+1):M:numel(q);
    phaseCoeffs = q(idx);
    j = 0:numel(phaseCoeffs)-1;
    nz = phaseCoeffs~=0;
    Gsp = sparse_terms(j(nz)*M,phaseCoeffs(nz));
    Fsp = sparse_terms(k,1);

    if isempty(Gsp.exponents)
        continue;
    end

    % DBB is checked from sparse supports before BMD multiplication.
    dbbOK = (max(double(Fsp.exponents)) < M) && ...
        all(mod(double(Gsp.exponents),M)==0) && numel(Gsp.exponents)>0;

    [mgr,refs] = build_bmd_pair_from_sparse_v10(Fsp,Gsp);
    nF = mgr.reachableNodeCount(refs(1,:));
    nG = mgr.reachableNodeCount(refs(2,:));
    st0 = mgr.stats();
    prod = mgr.multiply(refs(1,:),refs(2,:));
    st1 = mgr.stats(prod);
    nNew = st1.nodes_created - st0.nodes_created;
    nProd = mgr.reachableNodeCount(prod);

    expected = sparse_terms(double(Gsp.exponents)+k,Gsp.coefficients);
    got = sparse_terms_multiply(Fsp,Gsp);
    exactSparseOK = sparse_terms_same(got,expected);
    denseBMD = mgr.toDense(prod,1e5);
    denseExpected = sparse_terms_to_dense(expected);
    exactBMDOK = isequal(denseBMD,denseExpected);

    predNewOK = dbbOK && (nNew==nF);
    predSizeOK = dbbOK && (nProd==nF+nG);
    nontrivial = numel(Fsp.exponents)>=2 && nF>=2;

    branchRows = [branchRows; table(k,numel(Gsp.exponents),dbbOK,nontrivial,nF,nG,nNew,nProd, ...
        predNewOK,predSizeOK,exactSparseOK,exactBMDOK, ...
        'VariableNames',{'phase','phase_terms','dbb_precheck','nontrivial_low_factor', ...
        'V_F','V_G','N_new','V_product','prediction_new_ok','prediction_size_ok', ...
        'exact_sparse_ok','exact_bmd_ok'})]; %#ok<AGROW>
end
writetable(branchRows,fullfile(outDir,'polyphase_branch_results.csv'));

%% C. Summary with a deliberately conservative publication gate
strong = any(searchRows.dbb_factorization & searchRows.nontrivial_low_factor & ...
    searchRows.prediction_new_ok & searchRows.prediction_size_ok & searchRows.exact_product_ok);
branchExact = all(branchRows.dbb_precheck & branchRows.prediction_new_ok & ...
    branchRows.prediction_size_ok & branchRows.exact_sparse_ok & branchRows.exact_bmd_ok);

fid=fopen(fullfile(outDir,'summary.txt'),'w');
fprintf(fid,'Polyphase FIR DBB application validation\n');
fprintf(fid,'Source: MathWorks FIR Decimation for FPGA example\n');
fprintf(fid,'Design: firpm order 30, decimation factor 8, Q15 integer quantization\n');
fprintf(fid,'Whole-filter nontrivial DBB found: %d\n',strong);
fprintf(fid,'All natural polyphase branches satisfy exact DBB predictions: %d\n',branchExact);
if strong
    fprintf(fid,'Publication interpretation: STRONG external application candidate.\n');
elseif branchExact
    fprintf(fid,'Publication interpretation: STRUCTURAL SANITY CHECK ONLY; natural DBB cases are monomial low factors.\n');
else
    fprintf(fid,'Publication interpretation: DO NOT CLAIM application validation; inspect failures.\n');
end
fclose(fid);

disp(searchRows);
disp(branchRows);
fprintf('\nStrong external application candidate: %d\n',strong);
fprintf('All natural branches exact: %d\n',branchExact);

%% Local helper: exact integer rank-one factorization by residue modulo stride
function [ok,Fsp,Gsp,why] = exact_dbb_factorization(P,stride)
    maxE = double(max(P.exponents));
    a = zeros(1,maxE+1);
    a(double(P.exponents)+1) = P.coefficients;
    cols = ceil(numel(a)/stride);
    A = zeros(stride,cols);
    for e=0:maxE
        A(mod(e,stride)+1,floor(e/stride)+1)=a(e+1);
    end
    activeRows=find(any(A~=0,2));
    activeCols=find(any(A~=0,1));
    if numel(activeRows)<2 || numel(activeCols)<2
        ok=false; Fsp=sparse_terms([],[]); Gsp=sparse_terms([],[]);
        why='degenerate matrix; not a nontrivial two-sided factorization'; return;
    end
    r0=activeRows(1); c0=activeCols(1); pivot=A(r0,c0);
    % Exact rank-one test without floating division: A(r,c)*pivot=A(r,c0)*A(r0,c).
    for r=activeRows
        for c=activeCols
            if A(r,c)*pivot ~= A(r,c0)*A(r0,c)
                ok=false; Fsp=sparse_terms([],[]); Gsp=sparse_terms([],[]);
                why='coefficient matrix is not rank one'; return;
            end
        end
    end
    % Need an exact integer factorization compatible with BMD integer weights.
    row = A(:,c0).';
    g = 0;
    for v=row(row~=0), g=gcd(g,abs(round(v))); end
    if g==0
        ok=false; Fsp=sparse_terms([],[]); Gsp=sparse_terms([],[]);
        why='zero pivot column'; return;
    end
    f = row/g;
    nzf=find(f~=0,1);
    denom=f(nzf);
    h=zeros(1,cols);
    for c=activeCols
        val=A(nzf,c)/denom;
        if val~=round(val)
            ok=false; Fsp=sparse_terms([],[]); Gsp=sparse_terms([],[]);
            why='rank one over rationals but not exact integer factorization'; return;
        end
        h(c)=round(val);
    end
    fr=find(f~=0)-1; fc=f(f~=0);
    hc=find(h~=0)-1; hv=h(h~=0);
    Fsp=sparse_terms(fr,fc);
    Gsp=sparse_terms(hc*stride,hv);
    ok=true; why='exact integer DBB factorization';
end
