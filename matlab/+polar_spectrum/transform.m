function v=transform(v,transpose)
% Row-vector polar transform F^(tensor n), without bit reversal.
if nargin<2,transpose=false;end
N=numel(v);v=logical(v);step=1;
while step<N
    for start=1:2*step:N
        a=start:start+step-1;b=a+step;
        if transpose,v(b)=xor(v(b),v(a));else,v(a)=xor(v(a),v(b));end
    end
    step=2*step;
end
end
