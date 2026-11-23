classdef BigNat
% Exact nonnegative integers, little-endian base-10^7 limbs in doubles.
% Each multiplication/reduction stays below 2^53. No toolbox required.
methods (Static)
    function c=add(a,b)
        L=max(numel(a),numel(b));a(end+1:L)=0;b(end+1:L)=0;
        c=zeros(1,L);carry=0;
        for k=1:L
            v=a(k)+b(k)+carry;carry=floor(v/1e7);c(k)=v-carry*1e7;
        end
        if carry, c(end+1)=carry;end
    end
    function c=mul(a,b)
        if isempty(a)||isempty(b),c=[];return;end
        c=zeros(1,numel(a)+numel(b));
        for i=1:numel(a)
            carry=0;
            for j=1:numel(b)
                v=a(i)*b(j)+c(i+j-1)+carry;
                carry=floor(v/1e7);c(i+j-1)=v-carry*1e7;
            end
            c(i+numel(b))=carry;
        end
        while ~isempty(c)&&c(end)==0,c(end)=[];end
    end
    function a=pow2(k)
        a=1;for i=1:k,a=polar_spectrum.BigNat.add(a,a);end
    end
    function s=text(a)
        if isempty(a),s='0';return;end
        s=sprintf('%.0f',a(end));
        for k=numel(a)-1:-1:1,s=[s sprintf('%07.0f',a(k))];end %#ok<AGROW>
    end
end
end
