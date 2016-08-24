%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

clear all
clc

syms F u1x u1y u2x u2y f1x f1y f2x f2y v1x v1y v2x v2y

syms e1x e1y e2x e2y

syms traceE

syms mu lambda

F = [f1x f1y; f2x f2y];

E = [e1x e1y; e2x e2y];

I = eye(2,2);

GradU = [u1x u1y; u2x u2y ];

GradV = [v1x v1y; v2x v2y ];

FT = F.';

dE = 0.5 * (GradU.' * F + F.' * GradU);

% All
DF =    GradU * ( 2 * mu * E + lambda * traceE * I) ...
        + F * (2 * mu * dE + lambda * trace(dE) * I);

integral = sum(sum(DF.*GradV));


pos2der = {'x', 'y'};

for i = 1 : 2 % test
    for j = 1 : 2 % trial
        A{i,j} = integral;
        
        for ii = 1 : 2
            if ii ~= i
                for d = 1 : 2
                    variable = sprintf('v%d%s',ii,pos2der{d});
                    eval( ['A{i,j} = subs(A{i,j},', variable,', 0);'] );
                end
            end
        end
        
        for jj = 1 : 2
            if jj ~= j
                for d = 1 : 2
                    variable = sprintf('u%d%s',jj,pos2der{d});
                    eval( ['A{i,j} = subs(A{i,j},', variable,', 0);'] );
                end
            end
        end
        
    end
end

TableSubs = {'traceE', 'traceE[q]';...
             'f1x', 'F[q][0][0]';...
             'f1y', 'F[q][0][1]';...
             'f2x', 'F[q][1][0]';...
             'f2y', 'F[q][1][1]';...
             'e1x', 'E[q][0][0]';...
             'e1y', 'E[q][0][1]';...
             'e2x', 'E[q][1][0]';...
             'e2y', 'E[q][1][1]';...
             'u1x', 'gradphi[q][0][b]';...
             'u1y', 'gradphi[q][1][b]';...
             'u2x', 'gradphi[q][0][b]';...
             'u2y', 'gradphi[q][1][b]';...
             'v1x', 'gradphi[q][0][a]';...
             'v1y', 'gradphi[q][1][a]';...
             'v2x', 'gradphi[q][0][a]';...
             'v2y', 'gradphi[q][1][a]';...
             '/3', '/3.0';...
             '/2', '/2.0';...
             '/9', '/9.0';...
             '2*', '2.0*'};

for i = 1 : 2 % test
    for j = 1 : 2 % trial
        
        B{i,j} = A{i,j};
        fprintf('\n',i,j);%** Block %d%d ** \n
        
        
        fprintf('\naloc[a][%d][b][%d] += ( ',i-1,j-1);
        sA = size(B{i,j});
        for row=1:sA(1)
            for col = 1:sA(2)
                
                str2 = char(B{i,j}(row,col));
                
                for kk = 1 : size(TableSubs, 1)
                   
                    str2 = strrep(str2, TableSubs{kk,1}, TableSubs{kk,2});
                    
                end
                
                fprintf('%s ) * w[q]',str2);
            end
            %    fprintf(fid,';...\n');
        end
        fprintf(';');
        %fprintf('\niii = iii + 1;')
        
    end
end
fprintf('\n\n');
