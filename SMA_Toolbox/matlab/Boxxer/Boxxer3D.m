classdef Boxxer3D < Boxxer
    methods
        function obj = Boxxer3D( varargin )
            obj=obj@Boxxer(@Boxxer3D_Iface);
            if nargin>0
                obj.initialize(varargin{:});
            end
        end

%         function success=initialize(obj, imsize, sigma)
%             success=initialize@Boxxer(obj, imsize, sigma);
%         end

        function checkMaxima(obj, image, maxima, max_vals)
            Nmaxima=length(max_vals);
            for n=1:Nmaxima
                x=maxima(1,n);
                y=maxima(2,n);
                z=maxima(3,n);
                t=maxima(5,n);
                val=image(x,y,z,t);
                if val~=max_vals(n)
                    fprintf('(%i,%i,%i,%i):%.9f~=%.9f\n',x,y,z,t,val,max_vals(n));
                end
            end
        end
        
    end
end %classdef
