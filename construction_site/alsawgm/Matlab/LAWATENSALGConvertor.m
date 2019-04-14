classdef LAWATENSALGConvertor
    %
    % This file is part of the Matlab Toolbox TENSALG for Tensor Algebra,
    % developed under the BSD Licence.
    % See the LICENSE file for conditions.
    
    properties (Constant)
        storage = 'column'
    end
    
    methods (Static)
        
        function A = lawa2tensalgLaplaceOperator(Alawa)
            d = length(Alawa);
            n = cellfun(@(x) size(x,1),Alawa);
            As = cell(d,1);
            for mu = 1:d
                As{mu} = cell(d,1);
                for k = 1:d
                    if k == mu
                        As{mu}{k} = Alawa{mu};
                    else
                        As{mu}{k} = speye(n(mu));
                    end
                end
            end
            As = TSpaceOperators(As);
            Ac = DiagonalTensor(ones(d,1),d);
            A = TuckerLikeTensor(Ac,As);
        end
        
        function b = lawa2tensalgCanonicalTensor(blawa)
            d = length(blawa);
            n = cellfun(@(x) size(x,1),blawa);
            b =  TuckerLikeTensor(DiagonalTensor(1,d),TSpaceVectors(blawa));
            
        end
        
        function x = lawa2tensalgHTCores2TT(xlawa,r)
            d = length(xlawa)+1;
            n = cellfun(@(x) size(x,2),xlawa);
            
            x = cell(d,1);
            for mu=1:d-1
                switch LAWATENSALGConvertor.storage
                    case 'column'
                        x{mu+1} = FullTensor(xlawa{mu},3,[size(xlawa{mu},2),size(xlawa{mu},1)/r(mu),r(mu)]);
                    case 'row'
                        xlawa{mu}=xlawa{mu}';
                        x{mu+1} = FullTensor(xlawa{mu},[size(xlawa{mu},1),size(xlawa{mu},2)/r(mu),r(mu)]);
                        x{mu+1} = permute(x{mu+1},[2,1,3]);
                end
            end
            x{1} = FullTensor(eye(n(1)),3,[1,n(1),n(1)]);
            x = TTTensor(x);
            
        end
        
        function [xlawa,r] = tensalg2lawaTT2HTCores(x)
            
            d = x.order;
            x.cores{2} = squeeze(timesTensor(x.cores{1},x.cores{2},3,1),1);
            xlawa = x.cores(2:end);
            r = zeros(1,d-1);
            for mu=1:d-1
                r(mu) = size(xlawa{mu},3);
                switch LAWATENSALGConvertor.storage
                    case 'column'
                        xlawa{mu}=reshape(double(xlawa{mu}),xlawa{mu}.sz(2)*xlawa{mu}.sz(3), xlawa{mu}.sz(1));
                    case 'row'
                        xlawa{mu} = permute(xlawa{mu},[1,3,2]);
                        xlawa{mu} = reshape(double(xlawa{mu}),xlawa{mu}.sz(1)*xlawa{mu}.sz(3),xlawa{mu}.sz(2));
                end
            end
            
        end
        
        
        
    end
    
    
end
