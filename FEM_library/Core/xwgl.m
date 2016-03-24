function [x,w,flag] = xwgl(np,a,b)
%XWGL      Weights and nodes of the Gauss-Legendre
%           quadrature formula.
%
%    [X,W]=XWGL(NP) return the NP weighs and the nodes
%    of the corresponding Gauss-Legendre quadrature
%    formula in the reference interval (-1,1).
%
%    [X,W]=XWGL(NP,A,B) return the NP weighs and the nodes
%    of the corresponding Gauss-Legendre quadrature
%    formula in the reference interval (a,b).
%
%    Author: F. Saleri

if nargout == 3
    flag = 0;
end

n=np-1;
if np<=1
    x=0;w=2;
    return
end
x=zeros(np,1);
w=zeros(np,1);
jac=zeros(np);
k=[1:n];
v=(k)./(sqrt(4*(k.^2)-1));
jac=jac+diag(v,1)+diag(v,-1);
[w1,x]=eig(jac);
norm2=sqrt(diag(w1'*w1));    % normalizzazione per pesi
w1=(2*w1(1,:)'.^2)./norm2;   % 2=mu0; ' per aver w colonna
x=diag(x);		     % x colonna
%
% Ordino i nodi ed i rispettivi pesi
%
[x,ip]=sort(x);
for i=1:np
    w(i)=w1(ip(i));
end
%
% mappatura su (a,b)
%
if nargin == 3
    bma=(b-a)*.5;
    bpa=(b+a)*.5;
    x=bma*x+bpa;
    w=w*bma;
end

x = x';
w = w';

return