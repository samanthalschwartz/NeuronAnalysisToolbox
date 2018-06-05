

classdef DimerMCMC

    properties
        obsArr;  %cell array 1xNSeq.  Each cell is: t ax ay bx by SEax SEay SEbx SEby [nx9]
        trueState; %cell array 1xNSeq.  Each cell is Nx1 boolean 0=free, 1=bound
        NSeq;
        NObs; % array NSeq x 1
    end
    
    methods
        function obj = DimerMCMC(ObsData, TrueState)
            obj.obsArr = makecell(ObsData);
            obj.NSeq = numel(obj.obsArr);
            obj.NObs = cellfun(@(x) size(x,1),obj.obsArr);
            assert(all(cellfun(@(x) size(x,2),obj.obsArr)==9));
            if nargin>1
                obj.trueState = TrueState;
                assert(numel(obj.trueState)==obj.NSeq);
            end
        end
    end
    
end
