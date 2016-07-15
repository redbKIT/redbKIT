classdef HDF5_DenseMultiCVector < handle

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>    
    
    properties (GetAccess = public, SetAccess = protected)
        filename;
        dataset;
        rowDim;
        vecNumber;
        ComplexData;
    end
    
    methods
        %% Constructor
        function obj = HDF5_DenseMultiCVector(filename, dataset, n, isComplex)
            
            if nargin < 4
                isComplex = false;
            end
            obj.ComplexData = isComplex;
            obj.dataset   = dataset;
            obj.filename  = filename;
            
            % If file does not exists, create the file and the dataset
            if ~(exist(filename, 'file') == 2)
                
                if obj.ComplexData
                    h5create(filename, ['/',dataset,'/values_Re'],      [n Inf]  , 'Deflate', 0, 'ChunkSize', [n 1]);
                    h5create(filename, ['/',dataset,'/values_Im'],      [n Inf]  , 'Deflate', 0, 'ChunkSize', [n 1]);
                else
                    h5create(filename, ['/',dataset,'/values'],      [n Inf]  , 'Deflate', 0, 'ChunkSize', [n 1]);
                end
                h5create(filename, ['/',dataset,'/NumRows'],     [1 1]    );
                h5create(filename, ['/',dataset,'/NumVectors'],  [1 1]    );
                h5write(filename,  ['/',dataset,'/NumRows'],     int32(n) );
                
                if nargin < 3 || isempty( n )
                    obj.rowDim    = 0;
                else
                    obj.rowDim    = n;
                end
                
                obj.vecNumber = 0;
                
            else
                % file exists, dataset doesn't
                try
                    if obj.ComplexData
                        h5create(filename, ['/',dataset,'/values_Re'],      [n Inf]  , 'Deflate', 0, 'ChunkSize', [n 1]);
                        h5create(filename, ['/',dataset,'/values_Im'],      [n Inf]  , 'Deflate', 0, 'ChunkSize', [n 1]);
                    else
                        h5create(filename, ['/',dataset,'/values'],      [n Inf]  , 'Deflate', 0, 'ChunkSize', [n 1]);
                    end
                    h5create(filename, ['/',dataset,'/NumRows'],     [1 1]    );
                    h5create(filename, ['/',dataset,'/NumVectors'],  [1 1]    );
                    h5write(filename,  ['/',dataset,'/NumRows'],     int32(n) );
                    
                    obj.rowDim    = n;
                    obj.vecNumber = 0;
                    
                % file and dataset already exist
                catch
                    obj.rowDim    = h5read(obj.filename, ['/',obj.dataset,'/NumRows/']);
                    obj.vecNumber = h5read(obj.filename, ['/',obj.dataset,'/NumVectors/']);
                    if nargin==3 && n~=obj.rowDim
                        error('row dimension mismatch')
                    end
                end
            end
            
        end
        
        %% Append Vector
        function obj = append(obj,V)
            
            if issparse(V)
                V = full(V);
            end
            [n, m] = size(V);
            
            if n~=obj.rowDim
                error('Row dimension mismatch')
            end
            
            start = [1 obj.vecNumber+1];
            count = [n m];            
            if obj.ComplexData
                h5write(obj.filename, ['/',obj.dataset,'/values_Re'], real(V), start, count);
                h5write(obj.filename, ['/',obj.dataset,'/values_Im'], imag(V), start, count);
            else
                h5write(obj.filename, ['/',obj.dataset,'/values'], V, start, count);
            end
            h5write(obj.filename, ['/',obj.dataset,'/NumVectors'], int32(obj.vecNumber+m));
            
            obj.vecNumber = obj.vecNumber + m;
        end
        
        %% HDF5 info dataset
        function disp_dataset(obj)
            h5disp(obj.filename, ['/',obj.dataset]);
        end
        
        %% HDF5 info file
        function disp_file(obj)
            h5disp(obj.filename);
        end
        
        %% Read HDF5 file
        function V = readValues(obj, i)
            
            if nargin < 2
                
                if obj.ComplexData
                    V_Re  = h5read(obj.filename, ['/',obj.dataset,'/values_Re/']);
                    V_Im  = h5read(obj.filename, ['/',obj.dataset,'/values_Im/']);
                    V     = V_Re + complex(0,1) * V_Im;
                else
                    V  = h5read(obj.filename, ['/',obj.dataset,'/values/']);
                end
                
            end
            
            if nargin == 2
                
                if obj.ComplexData
                    V_Re  = h5read(obj.filename, ['/',obj.dataset,'/values_Re/'], [1 i], [obj.rowDim 1]);
                    V_Im  = h5read(obj.filename, ['/',obj.dataset,'/values_Im/'], [1 i], [obj.rowDim 1]);
                    V     = V_Re + complex(0,1) * V_Im;
                else
                    V  = h5read(obj.filename, ['/',obj.dataset,'/values/'], [1 i], [obj.rowDim 1]);
                end
                
            end
        end
        
    end
end
