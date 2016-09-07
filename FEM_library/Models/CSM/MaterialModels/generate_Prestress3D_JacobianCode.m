%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

clear all
clc

syms U u1 u2 u3 F u1x u1y u1z u2x u2y u2z u3x u3y u3z s1x s1y s1z s2x s2y s2z s3x s3y s3z v1x v1y v1z v2x v2y v2z v3x v3y v3z


S_0 = [s1x s1y s1z; s2x s2y s2z; s3x s3y s3z];

GradU = [u1x u1y u1z; u2x u2y u2z; u3x u3y u3z];

GradV = [v1x v1y v1z; v2x v2y v2z; v3x v3y v3z];

% All
DF =    GradU * S_0;

integral = sum(sum(DF.*GradV));


pos2der = {'x', 'y', 'z'};

for i = 1 : 3 % test
    for j = 1 : 3 % trial
        A{i,j} = integral;
        
        for ii = 1 : 3
            if ii ~= i
                for d = 1 : 3
                    variable = sprintf('v%d%s',ii,pos2der{d});
                    eval( ['A{i,j} = subs(A{i,j},', variable,', 0);'] );
                end
            end
        end
        
        for jj = 1 : 3
            if jj ~= j
                for d = 1 : 3
                    variable = sprintf('u%d%s',jj,pos2der{d});
                    eval( ['A{i,j} = subs(A{i,j},', variable,', 0);'] );
                end
            end
        end
        
    end
end

TableSubs = {'u1x', 'gradphi[q][0][b]';...
             'u1y', 'gradphi[q][1][b]';...
             'u1z', 'gradphi[q][2][b]';...
             'u2x', 'gradphi[q][0][b]';...
             'u2y', 'gradphi[q][1][b]';...
             'u2z', 'gradphi[q][2][b]';...
             'u3x', 'gradphi[q][0][b]';...
             'u3y', 'gradphi[q][1][b]';...
             'u3z', 'gradphi[q][2][b]';...
             's1x', 'PRESTRESS_NP1(ie,q,0,0)';...
             's1y', 'PRESTRESS_NP1(ie,q,0,1)';...
             's1z', 'PRESTRESS_NP1(ie,q,0,2)';...
             's2x', 'PRESTRESS_NP1(ie,q,1,0)';...
             's2y', 'PRESTRESS_NP1(ie,q,1,1)';...
             's2z', 'PRESTRESS_NP1(ie,q,1,2)';...
             's3x', 'PRESTRESS_NP1(ie,q,2,0)';...
             's3y', 'PRESTRESS_NP1(ie,q,2,1)';...
             's3z', 'PRESTRESS_NP1(ie,q,2,2)';...
             'v1x', 'gradphi[q][0][a]';...
             'v1y', 'gradphi[q][1][a]';...
             'v1z', 'gradphi[q][2][a]';...
             'v2x', 'gradphi[q][0][a]';...
             'v2y', 'gradphi[q][1][a]';...
             'v2z', 'gradphi[q][2][a]';...
             'v3x', 'gradphi[q][0][a]';...
             'v3y', 'gradphi[q][1][a]';...
             'v3z', 'gradphi[q][2][a]';...
             '/3', '/3.0';...
             '/2', '/2.0';...
             '/9', '/9.0';...
             '2*', '2.0*'};

for i = 1 : 3 % test
    for j = 1 : 3 % trial
        
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
