function n = workspace_bytes_v11(x)
%WORKSPACE_BYTES_V11 Bytes reported by WHOS for a concrete MATLAB value.
tmp=x; %#ok<NASGU>
w=whos('tmp');
n=double(w.bytes);
end
