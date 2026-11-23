function test_minimum
% Verify automatic minimum layers against the independent full small spectra.
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'matlab'));
vectors=jsondecode(fileread(fullfile(root,'tests','oracle_vectors.json')));
comparisons=0;
for k=1:numel(vectors)
    v=vectors(k);if v.D~=v.n,continue;end
    H=false(numel(v.rows),v.n);
    for j=1:v.n,H(:,j)=bitget(uint32(v.rows(:)),j);end
    weights=v.weights(:);counts=string(v.counts(:));
    for policy={'auto','enumerate','dp'}
        r=polar_spectrum.algorithm1(H,'auto',v.P(:)'+1,v.S(:)'+1,policy{1});
        if numel(weights)==1
            assert(isempty(r.dmin)&&strcmp(r.A_dmin,'0')&&isempty(r.search_cutoffs));
        else
            assert(r.dmin==weights(2)&&string(r.A_dmin)==counts(2));
            assert(isequal(r.weights(:),weights(1:2)));
            assert(r.search_cutoffs(end)==r.dmin && all(r.search_cutoffs<=r.dmin));
        end
        comparisons=comparisons+1;
    end
end
H=logical([1 1 0 0;1 0 1 0;0 0 0 1]);
r=polar_spectrum.algorithm1(H,'auto');
assert(r.dmin==3 && strcmp(r.A_dmin,'1') && isequal(r.search_cutoffs,[1 2 3]));
try
    polar_spectrum.algorithm1(H,'auto',[],[],'auto',r.stats.work_units-1);
    error('Expected cumulative work limit.');
catch e
    assert(strcmp(e.identifier,'polar_spectrum:workLimit'));
end
H=logical([eye(15) zeros(15,1)]);r=polar_spectrum.algorithm1(H,'auto');
assert(r.dmin==16 && strcmp(r.A_dmin,'1') && isequal(r.search_cutoffs,16));
report=struct('matlab_version',version,'oracle_comparisons',comparisons,'cases',[]);
for N=[64 128 256]
    for kind={'polar','pac'}
        H=polar_spectrum.checks_for_profile(polar_spectrum.profile(N,kind{1}));
        for mode={'none','puncture','shorten','mixed'}
            [P,S]=polar_spectrum.rate_matching(N,mode{1});
            truth=polar_spectrum.algorithm1(H,8,P,S);
            r=polar_spectrum.algorithm1(H,'auto',P,S);
            assert(r.dmin==truth.weights(2)&&strcmp(r.A_dmin,truth.counts{2}));
            r.kind=kind{1};r.mode=mode{1};report.cases=[report.cases;r]; %#ok<AGROW>
        end
    end
end
fid=fopen(fullfile(root,'results','matlab_minimum_validation.json'),'w');
cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
fprintf('PASS: %d automatic-minimum oracle comparisons; 24 profile cases; cumulative budget; dmin>8.\n',comparisons);
end
