function results = run_structural_representation_v04(varargin)
%RUN_STRUCTURAL_REPRESENTATION_V04 Scaling with explicit structural sharing.
%
% Grid family:
%   P = (sum_{i=0}^{innerN} x^i) *
%       (sum_{j=0}^{blocks-1} x^(j*2^stridePower)), innerN < 2^stridePower.
%
% Sparse storage grows with (innerN+1)*blocks terms.  The ordered BMD can
% reuse one high-bit block tail across the low-bit payload structure.
p=inputParser;
addParameter(p,'Quick',true,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'SaveResults',true,@(x)islogical(x)||ismember(x,[0 1]));
parse(p,varargin{:}); quick=logical(p.Results.Quick);
if quick
    cases=[31 7 8; 31 7 32; 31 7 128; 63 8 16; 63 8 64; 63 8 256; 255 10 16; 255 10 64; 255 10 256];
else
    cases=[15 6 4; 31 7 8; 31 7 16; 31 7 32; 31 7 64; 31 7 128; ...
           63 8 8; 63 8 16; 63 8 32; 63 8 64; 63 8 128; 63 8 256; ...
           255 10 8; 255 10 16; 255 10 32; 255 10 64; 255 10 128; 255 10 256];
end
rows=repmat(emptyRow(),size(cases,1),1);
fprintf('\nBMD-MATLAB v0.4 structural representation\n');
fprintf('=========================================\n');
for ii=1:size(cases,1)
    innerN=cases(ii,1); stridePower=cases(ii,2); blocks=cases(ii,3);
    row=emptyRow(); row.inner_n=innerN; row.stride_power=stridePower; row.blocks=blocks;
    stride=2^stridePower; row.degree=(blocks-1)*stride+innerN; row.logical_terms=(innerN+1)*blocks;
    fprintf('  inner=%-4g blocks=%-4g stride=2^%-2g ... ',innerN,blocks,stridePower);
    try
        row.bmd_build_time_s=timeit(@() bmdGridSink(innerN,blocks,stridePower));
        mb=BMDManager(); rb=mb.geometricGrid(innerN,blocks,stridePower); sb=mb.stats(rb);
        row.bmd_nodes=sb.reachable_nodes; row.bmd_workspace_nodes=sb.total_internal_nodes;
        row.bmd_garbage_nodes=sb.garbage_internal_nodes;
        row.bmd_payload_bytes_lb=sb.reachable_payload_bytes_lower_bound;
        if sb.garbage_internal_nodes~=0, error('BMD:GridGarbage','Direct grid builder created garbage nodes.'); end

        row.dense_build_time_s=timeit(@() denseGridSink(innerN,blocks,stridePower));
        dense=build_dense_grid(innerN,blocks,stridePower); wd=whos('dense');
        row.dense_coefficients=numel(dense); row.dense_bytes=wd.bytes;

        row.sparse_build_time_s=timeit(@() sparseGridSink(innerN,blocks,stridePower));
        sparse=build_sparse_grid(innerN,blocks,stridePower);
        row.sparse_terms=numel(sparse.exponents); row.sparse_bytes=sparse_terms_bytes(sparse);

        row.dense_bytes_over_bmd=row.dense_bytes/max(1,row.bmd_payload_bytes_lb);
        row.sparse_bytes_over_bmd=row.sparse_bytes/max(1,row.bmd_payload_bytes_lb);
        row.sparse_terms_per_bmd_node=row.sparse_terms/max(1,row.bmd_nodes);
        row.dense_build_over_bmd=row.dense_build_time_s/max(realmin,row.bmd_build_time_s);
        row.sparse_build_over_bmd=row.sparse_build_time_s/max(realmin,row.bmd_build_time_s);
        xv=0.999;
        yb=mb.evaluate(rb,xv); yd=polyval(dense,xv); ys=sparse_terms_evaluate(sparse,xv);
        row.numeric_check=near(yb,yd) && near(yb,ys);
        if ~row.numeric_check, error('BMD:GridNumeric','Grid numeric validation failed.'); end
        row.status='OK';
    catch ME
        row.status=['ERROR:' ME.identifier];
    end
    rows(ii)=row; fprintf('%s\n',row.status);
end
results=struct2table(rows); disp(results);
if p.Results.SaveResults
    outDir=fullfile(fileparts(mfilename('fullpath')),'results'); if ~exist(outDir,'dir'), mkdir(outDir); end
    writetable(results,fullfile(outDir,'structural_representation_v04.csv'));
end
end

function y=bmdGridSink(innerN,blocks,stridePower)
persistent m r
m=BMDManager(); r=m.geometricGrid(innerN,blocks,stridePower); s=m.stats(r);
y=double(s.total_internal_nodes)+double(r(1))+double(r(2));
end
function y=denseGridSink(innerN,blocks,stridePower)
persistent p
p=build_dense_grid(innerN,blocks,stridePower); y=double(numel(p))+p(1)+p(end);
end
function y=sparseGridSink(innerN,blocks,stridePower)
persistent s
s=build_sparse_grid(innerN,blocks,stridePower); y=double(numel(s.exponents))+s.coefficients(1)+s.coefficients(end);
end
function tf=near(a,b), tf=abs(a-b)<=1e-9*max(1,max(abs(a),abs(b))); end
function r=emptyRow()
r=struct('inner_n',NaN,'stride_power',NaN,'blocks',NaN,'degree',NaN,'logical_terms',NaN,'status','', ...
    'bmd_build_time_s',NaN,'bmd_nodes',NaN,'bmd_workspace_nodes',NaN,'bmd_garbage_nodes',NaN,'bmd_payload_bytes_lb',NaN, ...
    'dense_build_time_s',NaN,'dense_coefficients',NaN,'dense_bytes',NaN, ...
    'sparse_build_time_s',NaN,'sparse_terms',NaN,'sparse_bytes',NaN, ...
    'dense_bytes_over_bmd',NaN,'sparse_bytes_over_bmd',NaN,'sparse_terms_per_bmd_node',NaN, ...
    'dense_build_over_bmd',NaN,'sparse_build_over_bmd',NaN,'numeric_check',false);
end
