
classdef Boxxer2D < Boxxer
    properties (Access=protected)
    end
    methods
        function obj = Boxxer2D( varargin )
            obj=obj@Boxxer(@Boxxer2D_Iface);
            if nargin>0
                obj.initialize(varargin{:});
            end
        end

        function success=initialize(obj, imsize, sigma)
            if length(imsize)~=2
                error('Boxxer2D:initialize','Expected image size to have dimension 2, got: "%i"',length(imsize));
            end
            success=initialize@Boxxer(obj, imsize, sigma);
        end

        function checkMaxima(obj, image, maxima, max_vals)
            Nmaxima=length(max_vals);
            for n=1:Nmaxima
                x=maxima(1,n);
                y=maxima(2,n);
                t=maxima(4,n);
                val=image(x,y,t);
                if val~=max_vals(n)
                    fprintf('(%i,%i,%i):%.9f~=%.9f\n',x,y,t,val,max_vals(n));
                end
            end
        end
    end
end %classdef
