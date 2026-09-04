%% CIC application validation for exact nontrivial DBB closure
% External to the frozen v1.1.1 60-case benchmark.
%
% Application model: one FIR stage of a power-of-two CIC decimator.
% Standard CIC transfer function:
%   H_CIC(z) = [sum_{k=0}^{R*M-1} z^(-k)]^N
% so the equivalent cascade contains N identical all-ones FIR stages.
%
% We use R=256, M=1 and validate one exact FIR stage S_256(x), x=z^-1.
% IMPORTANT: we do NOT claim deployed CIC hardware explicitly implements the
% factorization below. It is an exact algebraic factorization of the standard
% equivalent FIR stage, chosen at the natural radix-16 exponent boundary.

clear; clc;
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
frozenDir = fullfile(repoRoot,'artifact','v1.1.1');
addpath(frozenDir);
outDir = fullfile(fileparts(mfilename('fullpath')),'results');
if ~exist(outDir,'dir'), mkdir(outDir); end

%% Fixed application parameters
R = 256;
M = 1;
N = 3;              % representative CIC order; validation is for one FIR stage
s = 4;
stride = 2^s;       % 16
L = R*M;            % 256

% One standard equivalent CIC FIR stage:
% S_256(x) = 1 + x + ... + x^255.
P = sparse_terms(0:L-1,ones(1,L));

% Exact radix-boundary factorization:
% S_256(x) = S_16(x) * S_16(x^16).
Fsp = sparse_terms(0:stride-1,ones(1,stride));
Gsp = sparse_terms((0:(L/stride-1))*stride,ones(1,L/stride));

%% DBB precondition from sparse supports
lowDegree = max(double(Fsp.exponents));
highStrideOK = all(mod(double(Gsp.exponents),stride)==0);
dbbOK = lowDegree < stride && highStrideOK && numel(Gsp.exponents)>1;
nontrivial = numel(Fsp.exponents)>=2;

%% Exact sparse identity
sparseProd = sparse_terms_multiply(Fsp,Gsp);
exactSparseOK = sparse_terms_same(sparseProd,P);

%% Frozen BMD implementation: operand-only manager, then multiply
[mgr,refs] = build_bmd_pair_from_sparse_v10(Fsp,Gsp);
nF = mgr.reachableNodeCount(refs(1,:));
nG = mgr.reachableNodeCount(refs(2,:));
st0 = mgr.stats();
prod = mgr.multiply(refs(1,:),refs(2,:));
st1 = mgr.stats(prod);
nNew = st1.nodes_created - st0.nodes_created;
nProd = mgr.reachableNodeCount(prod);

predNewOK = dbbOK && (nNew == nF);
predSizeOK = dbbOK && (nProd == nF+nG);

%% Exact BMD product validation
bmdDense = mgr.toDense(prod,1e6);
expectedDense = ones(1,L);
exactBMDOK = isequal(bmdDense,expectedDense);

%% Result
strong = dbbOK && nontrivial && exactSparseOK && exactBMDOK && predNewOK && predSizeOK;
T = table(R,M,N,s,stride,L,numel(Fsp.exponents),numel(Gsp.exponents), ...
    lowDegree,dbbOK,nontrivial,nF,nG,nNew,nProd,predNewOK,predSizeOK, ...
    exactSparseOK,exactBMDOK,strong, ...
    'VariableNames',{'R','M','N','s','stride','stage_length','F_terms','G_terms', ...
    'deg_F','dbb_precheck','nontrivial_low_factor','V_F','V_G','N_new','V_product', ...
    'prediction_new_ok','prediction_size_ok','exact_sparse_ok','exact_bmd_ok', ...
    'strong_application_candidate'});
writetable(T,fullfile(outDir,'cic_dbb_result.csv'));

fid=fopen(fullfile(outDir,'summary.txt'),'w');
fprintf(fid,'CIC DBB application validation\n');
fprintf(fid,'Application parameters: R=%d, M=%d, N=%d; one equivalent FIR stage tested\n',R,M,N);
fprintf(fid,'Stage: S_%d(x) = 1 + x + ... + x^%d\n',L,L-1);
fprintf(fid,'Factorization: S_%d(x) = S_%d(x) * S_%d(x^%d)\n',L,stride,L/stride,stride);
fprintf(fid,'F terms: %d; G terms: %d; deg(F)=%d < %d\n',numel(Fsp.exponents),numel(Gsp.exponents),lowDegree,stride);
fprintf(fid,'DBB precheck: %d\n',dbbOK);
fprintf(fid,'Exact sparse product: %d\n',exactSparseOK);
fprintf(fid,'|V(F)|=%d; |V(G)|=%d; N_new=%d; |V(FG)|=%d\n',nF,nG,nNew,nProd);
fprintf(fid,'N_new = |V(F)|: %d\n',predNewOK);
fprintf(fid,'|V(FG)| = |V(F)| + |V(G)|: %d\n',predSizeOK);
fprintf(fid,'Exact BMD product: %d\n',exactBMDOK);
fprintf(fid,'Strong nontrivial application candidate: %d\n',strong);
fprintf(fid,'Interpretation caveat: exact factorization of the standard equivalent FIR stage; not a claim about the internal deployed CIC hardware architecture.\n');
fclose(fid);

disp(T);
fprintf('\nStrong nontrivial CIC application candidate: %d\n',strong);
