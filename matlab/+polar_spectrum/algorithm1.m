function result=algorithm1(H,D,P,S,policy,maxWork)
% Exact Algorithm 1 for H*u'=0, x=u*F^(tensor n).
% P,S use one-based output coordinates. Counts are decimal strings.
% D='auto' returns dmin and A_dmin, stopping at the first positive layer.
% Pure MATLAB; no Python, MEX, Java, or Symbolic Math Toolbox.
if nargin<3,P=[];end
if nargin<4,S=[];end
if nargin<5,policy='auto';end
if nargin<6,maxWork=1e7;end
N=size(H,2);
if ~ismatrix(H)||N<1||2^round(log2(N))~=N||any(H(:)~=0 & H(:)~=1)
    error('polar_spectrum:input','H must be a binary matrix with power-of-two column count.');
end
minimum=(ischar(D)||(isstring(D)&&isscalar(D))) && strcmp(D,'auto');
if ~minimum,validateattributes(D,{'numeric'},{'scalar','integer','nonnegative','finite'});end
validateattributes(maxWork,{'numeric'},{'scalar','integer','positive','finite','<=',flintmax});
if ~ismember(string(policy),["auto","enumerate","dp"]),error('polar_spectrum:input','Unknown policy.');end
for positions={P,S}
    v=positions{1};
    if ~isnumeric(v)||any(~isfinite(v))||any(v~=fix(v))||any(v<1 | v>N)||numel(unique(v))~=numel(v)
        error('polar_spectrum:input','Coordinates must be distinct integers in [1,N].');
    end
end
if ~isempty(intersect(P,S)),error('polar_spectrum:input','P and S overlap.');end
start=tic;H=logical(H);C=canonical(H,N);
originalDimension=N-nnz(any(C,2));
for j=S(:)'
    v=false(1,N);v(j)=true;C=insert(C,polar_spectrum.transform(v,true));
end
shortenedDimension=N-nnz(any(C,2));free=find(~any(C,2))';
words=false(numel(free),N);
for k=1:numel(free)
    v=false(1,N);v(free(k))=true;
    for j=find(any(C,2))'
        if parity(v & C(j,:)),v(j)=~v(j);end
    end
    words(k,:)=polar_spectrum.transform(v);
end
projected=false(N,N);selected=false(1,numel(free));
for k=numel(free):-1:1
    v=words(k,:);v([P(:);S(:)])=false;
    [projected,selected(k)]=insert(projected,v);
end
for i=free(~selected),C(i,i)=true;end
transmittedDimension=nnz(selected);
requestedD=D;transmittedLength=N-numel(P)-numel(S);
if minimum,D=0;else,D=min(D,transmittedLength);end
q=ones(1,N);q(P)=0;
if minimum,q(S)=transmittedLength+1;else,q(S)=D+1;end
preprocessing=toc(start);stats=struct('work_units',0);
contexts=containers.Map('KeyType','char','ValueType','any');
free=find(~any(C,2))';searchCutoffs=[];
spectrum=cell(1,D+1);spectrum{1}=1;
if minimum
    bounds=arrayfun(@(i)lower(q,i),free);first=transmittedLength+1;
    if ~isempty(bounds),first=max(1,min(bounds));end
    for threshold=first:transmittedLength
        use('distance_trials');searchCutoffs(end+1)=threshold; %#ok<AGROW>
        D=threshold;spectrum=cell(1,D+1);spectrum{1}=1;
        for i=free(bounds<=threshold),spectrum=polyadd(spectrum,count(0,q,i));end
        if any(~cellfun(@isempty,spectrum(2:end))),break;end
    end
    assert(isempty(free)||any(~cellfun(@isempty,spectrum(2:end))), ...
        'Nontrivial transmitted code has no nonzero word.');
else
    for i=free,spectrum=polyadd(spectrum,count(0,q,i));end
end
elapsed=toc(start);
assert(isequal(spectrum{1},1),'Projection-injective invariant failed.');
weights=find(~cellfun(@isempty,spectrum))-1;
counts=cellfun(@polar_spectrum.BigNat.text,spectrum(weights+1),'UniformOutput',false);
result=struct('weights',weights,'counts',{counts},'exact',true,'n',N, ...
    'cutoff',requestedD,'evaluated_cutoff',D,'backend','matlab','policy',char(policy), ...
    'original_dimension',originalDimension,'shortened_dimension',shortenedDimension, ...
    'transmitted_dimension',transmittedDimension,'kernel_dimension',shortenedDimension-transmittedDimension, ...
    'preprocessing_seconds',preprocessing,'counting_seconds',elapsed-preprocessing, ...
    'elapsed_seconds',elapsed,'stats',stats);
if minimum
    result.dmin=[];result.A_dmin='0';result.search_cutoffs=searchCutoffs;
    if numel(weights)>1,result.dmin=weights(2);result.A_dmin=counts{2};end
end

    function use(name,amount)
        if nargin<2,amount=1;end
        stats.work_units=stats.work_units+amount;
        if stats.work_units>maxWork
            error('polar_spectrum:workLimit','Algorithm 1 work limit exceeded; increase maxWork or lower D.');
        end
        if ~isfield(stats,name),stats.(name)=0;end
        stats.(name)=stats.(name)+amount;
    end
    function out=polyadd(a,b)
        out=a;
        for ii=find(~cellfun(@isempty,b)),out{ii}=polar_spectrum.BigNat.add(out{ii},b{ii});end
    end
    function out=multiply(a,b)
        out=cell(1,D+1);ia=find(~cellfun(@isempty,a));ib=find(~cellfun(@isempty,b));
        for ii=ia
            for jj=ib(ib<=D+2-ii)
                use('coefficient_products');
                out{ii+jj-1}=polar_spectrum.BigNat.add(out{ii+jj-1},polar_spectrum.BigNat.mul(a{ii},b{jj}));
            end
        end
    end
    function pair=wef(costs,prefix)
        use('wef_nodes');
        if numel(costs)==1
            a=cell(1,D+1);b=a;a{1}=1;
            if costs(1)<=D,b{costs(1)+1}=1;end
            pair={a,b};return
        end
        paired=prefix(1:2*floor(numel(prefix)/2));e=paired(1:2:end);o=paired(2:2:end);
        f=wef(costs(1:2:end),xor(e,o));g=wef(costs(2:2:end),o);
        if mod(numel(prefix),2)
            bit=prefix(end);pair={multiply(f{bit+1},g{1}),multiply(f{2-bit},g{2})};
        else
            pair={polyadd(multiply(f{1},g{1}),multiply(f{2},g{2})), ...
                  polyadd(multiply(f{1},g{2}),multiply(f{2},g{1}))};
        end
    end
    function R=rows(offset,sizeLocal)
        key=sprintf('%d:%d',offset,sizeLocal);
        if isKey(contexts,key),R=contexts(key);return;end
        m=sizeLocal/2;indices=offset+m+(1:m);indices=indices(any(C(indices,:),2));
        R=struct('right',false(numel(indices),m),'left',C(indices,offset+(1:m)));
        for kk=1:numel(indices)
            R.right(kk,:)=polar_spectrum.transform(C(indices(kk),offset+m+(1:m)),true);
        end
        contexts(key)=R;
    end
    function F=family(a,s,costs,R)
        use('affine_families');m=numel(s);alpha=costs(1:m);beta=costs(m+1:end);
        delta=alpha+beta;delta(s)=abs(alpha(s)-beta(s));
        F=struct('m',m,'base',sum(min(alpha(s),beta(s))),'eta',s & alpha>beta, ...
            'free',delta==0,'pivot',false(1,m),'rows',false(m,m),'rhs',false(1,m));
        F.E=D-F.base;F.inconsistent=F.E<0;res=false(m,m);resrhs=false(1,m);resused=false(1,m);
        for rr=1:size(R.right,1)
            row=R.right(rr,:);rhs=xor(parity(R.left(rr,:) & a),parity(row & F.eta));stored=false;
            while any(row & F.free)
                pp=find(row & F.free,1,'last');
                if ~F.pivot(pp),F.pivot(pp)=true;F.rows(pp,:)=row;F.rhs(pp)=rhs;stored=true;break;end
                row=xor(row,F.rows(pp,:));rhs=xor(rhs,F.rhs(pp));
            end
            if stored,continue;end
            while any(row)
                pp=find(row,1,'last');
                if ~resused(pp),resused(pp)=true;res(pp,:)=row;resrhs(pp)=rhs;stored=true;break;end
                row=xor(row,res(pp,:));rhs=xor(rhs,resrhs(pp));
            end
            if ~stored && rhs,F.inconsistent=true;end
        end
        F.rank=nnz(resused);F.target=resrhs(resused);F.factor=polar_spectrum.BigNat.pow2(nnz(F.free)-nnz(F.pivot));
        positive=find(delta>0);valuesPositive=delta(positive);
        [~,order]=sortrows([valuesPositive(:) positive(:)],[1 2]);F.indices=reshape(positive(order),1,[]);
        F.delta=delta(F.indices);F.columns=res(resused,F.indices)';
    end
    function patterns(F,emit)
        if F.inconsistent,return;end
        eligible=find(F.delta<=F.E);suffix=cell(1,numel(eligible)+1);suffix{end}=false(F.rank,F.rank);
        for kk=numel(eligible):-1:1,suffix{kk}=insert(suffix{kk+1},F.columns(eligible(kk),:));end
        deviations(F,eligible,suffix,1,F.target,false(1,F.m),0,emit);
    end
    function deviations(F,eligible,suffix,next,syndrome,error,cost,emit)
        use('deviation_nodes');
        if any(reduce(suffix{next},syndrome)),return;end
        if ~any(syndrome),emit(error,cost);end
        for kk=next:numel(eligible)
            item=eligible(kk);following=cost+F.delta(item);
            if following>F.E,break;end
            v=error;v(F.indices(item))=true;
            deviations(F,eligible,suffix,kk+1,xor(syndrome,F.columns(item,:)),v,following,emit);
        end
    end
    function out=histogram(F)
        out=cell(1,D+1);if F.inconsistent,return;end
        p=numel(F.delta);h=0;if p,h=min(p,floor(F.E/F.delta(1)));end
        qb=1;choose=1;for kk=1:h,choose=choose*(p-kk+1)/kk;qb=qb+choose;end
        dp=strcmp(policy,'dp') || (strcmp(policy,'auto') && (p+1)*(F.E+1)*2^F.rank<F.m^2*qb);
        if dp
            use('dp_families');states=containers.Map('KeyType','char','ValueType','any');
            states(statekey(0,false(1,F.rank)))=1;
            for kk=1:p
                keys0=keys(states);use('dp_states',numel(keys0));
                values0=values(states);following=containers.Map(keys0,values0);
                for jj=1:numel(keys0)
                    [cost,syndrome]=parsestate(keys0{jj});newcost=cost+F.delta(kk);
                    if newcost<=F.E
                        key=statekey(newcost,xor(syndrome,F.columns(kk,:)));
                        value=values0{jj};if isKey(following,key),value=polar_spectrum.BigNat.add(following(key),value);end
                        following(key)=value;
                    end
                end
                states=following;
            end
            keys0=keys(states);
            for kk=1:numel(keys0)
                [cost,syndrome]=parsestate(keys0{kk});
                if isequal(syndrome,F.target),out{F.base+cost+1}=polar_spectrum.BigNat.mul(F.factor,states(keys0{kk}));end
            end
        else
            use('enumerated_families');
            % Accumulate only this histogram, while Gen below remains streamed.
            collector=polar_spectrum.PolyCollector(D);
            patterns(F,@(~,cost)collector.add(F.base+cost,F.factor));out=collector.coefficients;
        end
    end
    function generate(offset,costs,first,emit)
        use('generation_calls');
        if lower(costs,first)>D,return;end
        if numel(costs)==1,use('generated_pairs');emit(true,true);return;end
        m=numel(costs)/2;
        if first>m
            generate(offset+m,combine(costs,true),first-m,@(b,t)repeatpair(b,t,emit));return
        end
        generate(offset,combine(costs,false),first,@(a,s)expand(a,s,costs,rows(offset,numel(costs)),emit));
    end
    function repeatpair(b,t,emit)
        use('generated_pairs');emit([false(size(b)) b],[t t]);
    end
    function expand(a,s,costs,R,emit)
        F=family(a,s,costs,R);independent=find(F.free & ~F.pivot);du=false(numel(independent),2*F.m);dx=du;
        for kk=1:numel(independent)
            v=false(1,F.m);v(independent(kk))=true;
            for pp=find(F.pivot)
                if parity(F.rows(pp,:) & F.free & v),v(pp)=~v(pp);end
            end
            du(kk,F.m+1:end)=polar_spectrum.transform(v);dx(kk,:)=[v v];
        end
        patterns(F,@(error,~)emitfamily(F,error,a,s,du,dx,emit));
    end
    function emitfamily(F,error,a,s,du,dx,emit)
        for pp=find(F.pivot)
            if parity(F.rows(pp,:) & error)~=F.rhs(pp),error(pp)=~error(pp);end
        end
        t=xor(error,F.eta);u=[a polar_spectrum.transform(t)];x=[xor(s,t) t];
        gray(u,x,size(du,1),du,dx,emit);
    end
    function [u,x]=gray(u,x,k,du,dx,emit)
        if k==0,use('generated_pairs');emit(u,x);return;end
        [u,x]=gray(u,x,k-1,du,dx,emit);u=xor(u,du(k,:));x=xor(x,dx(k,:));
        [u,x]=gray(u,x,k-1,du,dx,emit);
    end
    function out=count(offset,costs,first)
        use('count_calls');out=cell(1,D+1);
        if lower(costs,first)>D,return;end
        if ~any(any(C(offset+first+1:offset+numel(costs),:),2))
            use('wef_terminals');pair=wef(costs,false(1,first-1));out=pair{2};return
        end
        m=numel(costs)/2;
        if first>m,out=count(offset+m,combine(costs,true),first-m);return;end
        collector=polar_spectrum.PolyCollector(D);
        generate(offset,combine(costs,false),first,@(a,s)collector.update(histogram(family(a,s,costs,rows(offset,numel(costs))))));
        out=collector.coefficients;
    end
end

function C=canonical(H,N)
C=false(N,N);for k=1:size(H,1),C=insert(C,H(k,:));end
end
function [C,added]=insert(C,row)
added=false;
while any(row)
    j=find(row,1,'last');
    if ~C(j,j),C(j,:)=row;added=true;return;end
    row=xor(row,C(j,:));
end
end
function row=reduce(C,row)
for j=size(C,1):-1:1,if row(j)&&C(j,j),row=xor(row,C(j,:));end,end
end
function bit=parity(row)
bit=mod(nnz(row),2)~=0;
end
function costs=combine(q,repeat)
m=numel(q)/2;
if repeat,costs=q(1:m)+q(m+1:end);else,costs=min(q(1:m),q(m+1:end));end
end
function bound=lower(costs,first)
while numel(costs)>1
    m=numel(costs)/2;costs=combine(costs,first>m);first=mod(first-1,m)+1;
end
bound=costs(1);
end
function key=statekey(cost,syndrome)
key=[sprintf('%d:',cost) char('0'+syndrome)];
end
function [cost,syndrome]=parsestate(key)
separator=find(key==':',1);cost=str2double(key(1:separator-1));syndrome=key(separator+1:end)=='1';
end
