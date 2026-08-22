function write_run_metadata_v11(varargin)
%WRITE_RUN_METADATA_V11 Reproducibility metadata without machine identity.
p=inputParser; addParameter(p,'Trials',5); addParameter(p,'BuildTrials',3); addParameter(p,'DenseMaxCoefficients',250000); parse(p,varargin{:});
root=fileparts(mfilename('fullpath')); outDir=fullfile(root,'results'); if ~exist(outDir,'dir'), mkdir(outDir); end
fid=fopen(fullfile(outDir,'run_metadata_v11.txt'),'w'); if fid<0, error('BMD:V11Metadata','Cannot open metadata file.'); end
cl=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'BMD-MATLAB v1.1 publication-readiness validation\n');
fprintf(fid,'Generated: %s\n',datestr(now,'yyyy-mm-dd HH:MM:SS'));
fprintf(fid,'MATLAB version: %s\n',version);
fprintf(fid,'Computer: %s\n',computer);
fprintf(fid,'Frozen predictor: v0.9-static-20260820 (no threshold retuning).\n');
fprintf(fid,'Validation cases: 60 across 9 families; 54 application-formula cases + 6 controls.\n');
fprintf(fid,'Operation trials per case: %d\n',round(p.Results.Trials));
fprintf(fid,'BMD conversion trials per case: %d\n',round(p.Results.BuildTrials));
fprintf(fid,'Dense result coefficient cap: %d\n',round(p.Results.DenseMaxCoefficients));
fprintf(fid,'Sparse baselines: unique+accumarray and independent sort-reduce; best median is used for BMD-vs-sparse claims.\n');
fprintf(fid,'BMD memory metric: WHOS bytes of packed reachable numeric DAG only; live containers.Map/cache overhead excluded.\n');
fprintf(fid,'Dense baseline is diagnostic when result coefficient count is within the frozen cap.\n');
fprintf(fid,'Same package should be run unchanged on a second MATLAB environment before final submission.\n');
end
