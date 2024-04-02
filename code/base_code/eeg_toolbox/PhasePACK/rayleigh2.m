function [p,rBar]=rayleigh2(data,dim)
%function [p,rBar]=rayleigh2(data,dim)
%this is a vectorized version of rayleigh

if ~exist('dim','var'), dim=1; end

[data,perm,nshifts]=shiftdata(data,dim);

c=cos(data); s=sin(data);
n=size(c,1);
rBar=sqrt(sum(c,1).^2+sum(s,1).^2)./n;

z=n.*rBar.^2;
p = exp(-z) .* (1 + (2*z - z.^2) ./ (4.*n) - (24*z - 132*z.^2 + 76*z.^3 - 9*z.^4) ./ (288*n.^2));


p=unshiftdata(p,perm,nshifts);
if nargout==2
  rBar=unshiftdata(rBar,perm,nshifts);
end

