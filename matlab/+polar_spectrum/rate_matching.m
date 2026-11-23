function [P,S,ET,EP,ES] = rate_matching(N, mode)
% Fixed coordinate-selector matrices; x_tx=x*ET.', ES*x.'=0.
if nargin<2, mode='none'; end
P=[]; S=[];
switch N
case 64
switch char(mode)
case 'none'
P=[]; S=[];
case 'puncture'
P=[1 2 3 4]; S=[];
case 'shorten'
P=[]; S=[61 62 63 64];
case 'mixed'
P=[1 2]; S=[63 64];
otherwise, error('Unknown rate-matching mode.');
end
case 128
switch char(mode)
case 'none'
P=[]; S=[];
case 'puncture'
P=[1 2 3 4 5 6 7 8]; S=[];
case 'shorten'
P=[]; S=[121 122 123 124 125 126 127 128];
case 'mixed'
P=[1 2 3 4]; S=[125 126 127 128];
otherwise, error('Unknown rate-matching mode.');
end
case 256
switch char(mode)
case 'none'
P=[]; S=[];
case 'puncture'
P=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]; S=[];
case 'shorten'
P=[]; S=[241 242 243 244 245 246 247 248 249 250 251 252 253 254 255 256];
case 'mixed'
P=[1 2 3 4 5 6 7 8]; S=[249 250 251 252 253 254 255 256];
otherwise, error('Unknown rate-matching mode.');
end
otherwise, error('Supported lengths: 64,128,256.');
end
E=eye(N); EP=E(P,:); ES=E(S,:); ET=E(setdiff(1:N,[P S]),:);
end
