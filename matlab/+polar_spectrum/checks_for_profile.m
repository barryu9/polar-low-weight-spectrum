function H=checks_for_profile(cfg)
% Causal checks H*u'=0 from frozen source v and u=v*g.
N=cfg.n;g=cfg.taps;h=false(1,N);h(1)=true;
for i=2:N
    for lag=1:min(i-1,numel(g)-1)
        if g(lag+1),h(i)=xor(h(i),h(i-lag));end
    end
end
frozen=setdiff(1:N,cfg.information);H=false(numel(frozen),N);
for row=1:numel(frozen),i=frozen(row);H(row,1:i)=h(i:-1:1);end
end
