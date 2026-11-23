function result=run_spectrum(N,kind,mode,D,policy,maxWork)
% RUN_SPECTRUM(128,'pac','puncture',8) uses only the MATLAB implementation.
% RUN_SPECTRUM(128,'pac','puncture','auto') finds only dmin and A_dmin.
if nargin<1,N=64;end
if nargin<2,kind='polar';end
if nargin<3,mode='none';end
if nargin<4,D=8;end
if nargin<5,policy='auto';end
if nargin<6,maxWork=1e7;end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(fullfile(root,'matlab'));
cfg=polar_spectrum.profile(N,kind);
H=polar_spectrum.checks_for_profile(cfg);
[P,S,ET,EP,ES]=polar_spectrum.rate_matching(N,mode); %#ok<ASGLU>
result=polar_spectrum.algorithm1(H,D,P,S,policy,maxWork);
fprintf('%s mother (%d,%d), GA Eb/N0=2 dB, %s, MATLAB\n',upper(kind),N,N/2,mode);
if isfield(result,'dmin')
    fprintf('Transmitted dimension=%d; automatic dmin; exact=%d\n',result.transmitted_dimension,result.exact);
    if isempty(result.dmin),fprintf('dmin = undefined (zero-dimensional code)\n');
    else,fprintf('dmin = %d\n',result.dmin);end
    fprintf('A_dmin = %s\n',result.A_dmin);
else
    fprintf('Transmitted dimension=%d; D=%d; exact=%d\n',result.transmitted_dimension,D,result.exact);
    fprintf('weight   multiplicity\n');
    for k=1:numel(result.weights),fprintf('%6d   %s\n',result.weights(k),result.counts{k});end
end
fprintf('Algorithm time: %.6f s\n',result.elapsed_seconds);
end
