function [x,y] = mesh_square(bs,s)

nbs = 4;

if nargin==0
    x=nbs;
    return
end

mu1=1;


dl(:,1) = [0.000000; 1.000000; 1; 0];
dl(:,2) = [0.000000; 1.000000; 1; 0];
dl(:,3) = [0.000000; 1.000000; 1; 0];
dl(:,4) = [0.000000; 1.000000; 1; 0];


bs1=bs(:)';

if nargin==1
    x=dl(:,bs1);
    return
end

x=zeros(size(s));
y=zeros(size(s));
[m,n]=size(bs);
if m==1 && n==1,
    bs=bs*ones(size(s)); % expand bs
elseif m~=size(s,1) || n~=size(s,2),
    error('PDE:EllipsoidA2_pdedom:SizeBs', 'bs must be scalar or of same size as s.');
end

psym = [0, 0; mu1, 0; mu1, mu1; 0, mu1];


if ~isempty(s)
   % boundary segment 1
   ii=find(bs==1);
   if length(ii)
       t = s(ii);
       x(ii) = psym(1,1) + t*(psym(2,1)-psym(1,1));
       y(ii) = psym(1,2) + t*(psym(2,2)-psym(1,2));
   end
   
   % boundary segment 2
   ii=find(bs==2);
   if length(ii)
       t = s(ii);
       x(ii) = psym(2,1) + t*(psym(3,1)-psym(2,1));
       y(ii) = psym(2,2) + t*(psym(3,2)-psym(2,2));
   end
   
   % boundary segment 3
   ii=find(bs==3);
   if length(ii)
       t = s(ii);
       x(ii) = psym(3,1) + t*(psym(4,1)-psym(3,1));
       y(ii) = psym(3,2) + t*(psym(4,2)-psym(3,2));
   end

   % boundary segment 4
   ii=find(bs==4);
   if length(ii)
       t = s(ii);
       x(ii) = psym(4,1) + t*(psym(1,1)-psym(4,1));
       y(ii) = psym(4,2) + t*(psym(1,2)-psym(4,2));
   end
   
   
end

