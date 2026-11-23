function test_matlab
% Independent oracle checks, followed by all 24 fixed-profile examples.
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'matlab'));
cases=jsondecode(fileread(fullfile(root,'tests','oracle_vectors.json')));
policies={'auto','enumerate','dp'};comparisons=0;
for k=1:numel(cases)
    v=cases(k);H=false(numel(v.rows),v.n);
    if ~isempty(v.rows)
        for j=1:v.n,H(:,j)=bitget(uint32(v.rows(:)),j);end
    end
    for p=1:numel(policies)
        r=polar_spectrum.algorithm1(H,v.D,v.P(:)'+1,v.S(:)'+1,policies{p});
        assert(isequal(r.weights(:),v.weights(:)),sprintf('Weight mismatch in case %d/%s',k,policies{p}));
        assert(isequal(string(r.counts(:)),string(v.counts(:))),sprintf('Count mismatch in case %d/%s',k,policies{p}));
        comparisons=comparisons+1;
    end
end
try
    polar_spectrum.algorithm1(false(0,8),8,[],[],'auto',1);
    error('Expected work limit.');
catch e
    assert(strcmp(e.identifier,'polar_spectrum:workLimit'));
end
report=struct('matlab_version',version,'oracle_comparisons',comparisons,'cases',[]);
fixtures=jsondecode(fileread(fullfile(root,'profiles','ga_2db.json')));
for N=[64 128 256]
    for kind={'polar','pac'}
        cfg=polar_spectrum.profile(N,kind{1});H=polar_spectrum.checks_for_profile(cfg);
        fixture=fixtures.(sprintf('%s_%d_%d',kind{1},N,N/2));
        assert(isequal(cfg.information_zero_based(:),fixture.information(:)));
        assert(isequal(cfg.taps(:),fixture.taps(:)));
        for mode={'none','puncture','shorten','mixed'}
            [P,S,ET,EP,ES]=polar_spectrum.rate_matching(N,mode{1});
            assert(size(ET,1)+size(EP,1)+size(ES,1)==N);
            for which={'punctured','shortened'}
                selector=EP;if strcmp(which{1},'shortened'),selector=ES;end
                if ~isempty(selector)
                    literal=readmatrix(fullfile(root,'profiles','matrices',sprintf('%d_%s_%s.csv',N,mode{1},which{1})));
                    assert(isequal(selector,literal));
                end
            end
            r=polar_spectrum.algorithm1(H,8,P,S);
            r.kind=kind{1};r.mode=mode{1};
            report.cases=[report.cases;r]; %#ok<AGROW>
        end
    end
end
fid=fopen(fullfile(root,'results','matlab_validation.json'),'w');
cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
fprintf('PASS: %d oracle/policy comparisons; 24 profile cases saved.\n',comparisons);
end
