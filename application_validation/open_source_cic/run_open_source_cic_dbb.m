%% Open-source FPGA CIC application validation
% External application-derived structural validation.
% Upstream hardware-test configuration:
% davemuscle/sigma_delta_converters @
% 198bdecf66dd147b26ec9f4196e8bb03c9abfb53
% rtl/hw_test/top.sv: OSR=128, CIC=2.

clear; clc;

thisFile = mfilename('fullpath');
thisDir = fileparts(thisFile);
repoRoot = fileparts(fileparts(thisDir));
artifactDir = fullfile(repoRoot, 'artifact', 'v1.1.1');
addpath(artifactDir);

outDir = fullfile(thisDir, 'results');
if ~exist(outDir, 'dir'), mkdir(outDir); end

upstream_repo = "davemuscle/sigma_delta_converters";
upstream_commit = "198bdecf66dd147b26ec9f4196e8bb03c9abfb53";
upstream_file = "rtl/hw_test/top.sv";
R = 128;                 % upstream OSR
N = 2;                   % upstream CIC stages
M = 1;                   % standard comb differential delay for this RTL
s = 3;
stride = 2^s;

% Exact FIR-equivalent stage S_128(x).
S = sparse_terms(uint64(0:R-1), ones(1,R));

% Exact DBB factorization S_128(x) = S_8(x) S_16(x^8).
F = sparse_terms(uint64(0:stride-1), ones(1,stride));
G = sparse_terms(uint64(0:stride:R-1), ones(1,R/stride));

% DBB precondition.
degF = double(max(F.exponents));
gSupp = double(G.exponents);
gNonzero = gSupp(gSupp ~= 0);
if isempty(gNonzero)
    support_gcd = 0;
else
    support_gcd = gNonzero(1);
    for k = 2:numel(gNonzero)
        support_gcd = gcd(support_gcd, gNonzero(k));
    end
end
dbb_precheck = (degF < stride) && all(mod(gSupp, stride) == 0);

% Exact sparse reference multiplication.
prodExp = double(F.exponents(:)) + double(G.exponents(:)).';
prodCoef = double(F.coefficients(:)) * double(G.coefficients(:)).';
prodExp = prodExp(:);
prodCoef = prodCoef(:);
[uExp,~,ic] = unique(prodExp);
uCoef = accumarray(ic, prodCoef);
keep = (uCoef ~= 0);
P_sparse = sparse_terms(uint64(uExp(keep).'), uCoef(keep).');
exact_sparse_ok = isequal(P_sparse.exponents, S.exponents) && ...
    isequal(P_sparse.coefficients, S.coefficients);

% Build operands in the frozen implementation, compact to operand-only manager,
% then measure multiplication-induced closure. The frozen helper returns a
% manager and a 2-row refs array, not three separate outputs.
[mgr, refs] = build_bmd_pair_from_sparse_v10(F, G);
p = refs(1,:);
q = refs(2,:);
V_F = mgr.reachableNodeCount(p);
V_G = mgr.reachableNodeCount(q);
V_before = mgr.reachableNodeCount([p;q]);
next_before = mgr.stats().nextNodeId;
r = mgr.multiply(p, q);
next_after = mgr.stats().nextNodeId;
N_new = double(next_after - next_before);
V_product = mgr.reachableNodeCount(r);

% Exact product check against S_128.
S_dense = zeros(1,R);
S_dense(double(S.exponents)+1) = double(S.coefficients);
P_dense = mgr.toDense(r);
P_dense = P_dense(:).';
if numel(P_dense) < R, P_dense(end+1:R) = 0; end
exact_bmd_ok = isequal(P_dense(1:R), S_dense) && ...
    all(P_dense(R+1:end) == 0);

prediction_new_ok = (N_new == V_F);
prediction_size_ok = (V_product == V_F + V_G);
nontrivial_low_factor = (numel(F.exponents) >= 2) && (V_F >= 2);
strong_application_candidate = dbb_precheck && nontrivial_low_factor && ...
    exact_sparse_ok && exact_bmd_ok && prediction_new_ok && prediction_size_ok;

T = table(upstream_repo, upstream_commit, upstream_file, R, M, N, s, stride, ...
    numel(F.exponents), numel(G.exponents), degF, support_gcd, dbb_precheck, ...
    nontrivial_low_factor, V_F, V_G, V_before, N_new, V_product, ...
    prediction_new_ok, prediction_size_ok, exact_sparse_ok, exact_bmd_ok, ...
    strong_application_candidate, ...
    'VariableNames', {'upstream_repo','upstream_commit','upstream_file','R','M','N', ...
    's','stride','F_terms','G_terms','deg_F','support_gcd_G','dbb_precheck', ...
    'nontrivial_low_factor','V_F','V_G','V_before','N_new','V_product', ...
    'prediction_new_ok','prediction_size_ok','exact_sparse_ok','exact_bmd_ok', ...
    'strong_application_candidate'});

disp(T);
fprintf('Strong open-source FPGA CIC application candidate: %d\n', strong_application_candidate);

writetable(T, fullfile(outDir, 'open_source_cic_dbb_result.csv'));
fid = fopen(fullfile(outDir, 'summary.txt'), 'w');
fprintf(fid, 'Open-source FPGA CIC DBB application validation\n');
fprintf(fid, 'Upstream: %s @ %s\n', upstream_repo, upstream_commit);
fprintf(fid, 'Hardware-test config: R/OSR=%d, M=%d, N/CIC_STAGES=%d\n', R, M, N);
fprintf(fid, 'Exact DBB factorization: S_128(x)=S_8(x) S_16(x^8), s=%d, stride=%d\n', s, stride);
fprintf(fid, '|V(F)|=%d, |V(G)|=%d, N_new=%d, |V(FG)|=%d\n', V_F, V_G, N_new, V_product);
fprintf(fid, 'prediction_new_ok=%d\n', prediction_new_ok);
fprintf(fid, 'prediction_size_ok=%d\n', prediction_size_ok);
fprintf(fid, 'exact_sparse_ok=%d\n', exact_sparse_ok);
fprintf(fid, 'exact_bmd_ok=%d\n', exact_bmd_ok);
fprintf(fid, 'strong_application_candidate=%d\n', strong_application_candidate);
fprintf(fid, ['Interpretation: exact validation on an FIR-equivalent stage derived from ', ...
    'a public open-source FPGA hardware-test CIC configuration; no claim is made ', ...
    'that the RTL internally implements this factorization.\n']);
fclose(fid);
