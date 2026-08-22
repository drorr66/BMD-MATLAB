function write_run_metadata_v111(varargin)
%WRITE_RUN_METADATA_V111 Reproducibility metadata for validation-only patch.
p=inputParser; addParameter(p,'Trials',5); addParameter(p,'BuildTrials',3); addParameter(p,'DenseMaxCoefficients',250000); parse(p,varargin{:});
root=fileparts(mfilename('fullpath')); outDir=fullfile(root,'results'); if ~exist(outDir,'dir'), mkdir(outDir); end
fid=fopen(fullfile(outDir,'run_metadata_v111.txt'),'w'); if fid<0, error('BMD:V111Metadata','Cannot open metadata file.'); end
cl=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'BMD-MATLAB v1.1.1 publication validation-only patch\n');
fprintf(fid,'Generated: %s\n',datestr(now,'yyyy-mm-dd HH:MM:SS'));
fprintf(fid,'MATLAB version: %s\n',version);
fprintf(fid,'Computer: %s\n',computer);
fprintf(fid,'Frozen predictor: v0.9-static-20260820 (no threshold retuning).\n');
fprintf(fid,'Frozen workload: exact same 60 v1.1 cases across 9 families.\n');
fprintf(fid,'Operation trials per case: %d\n',round(p.Results.Trials));
fprintf(fid,'BMD conversion trials per case: %d\n',round(p.Results.BuildTrials));
fprintf(fid,'Dense result coefficient cap: %d\n',round(p.Results.DenseMaxCoefficients));
fprintf(fid,'Sparse baselines unchanged: unique+accumarray and sort-reduce; best median used for claims.\n');
fprintf(fid,'Timing methodology unchanged from v1.1; validation remains outside timed regions.\n');
fprintf(fid,'Validation patch: removed x=0.999 gate; dense-eligible BMD results use coefficient-for-coefficient comparison; dense-ineligible results use 3 deterministic exact-modular fingerprints.\n');
fprintf(fid,'Validation subchecks are emitted separately in publication_validation_results_v111.csv.\n');
fprintf(fid,'BMD memory metric unchanged: WHOS bytes of packed reachable numeric DAG only; live containers.Map/cache overhead excluded.\n');
fprintf(fid,'After a clean first-environment rerun, this exact ZIP must be run unchanged on a second MATLAB environment before final submission.\n');
end
