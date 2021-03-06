classdef Gauss2DMAP < MappelBase
    properties 
        Name='Gauss2DMAP';
        nParams=4;
        ParamNames={'x', 'y', 'I', 'bg'};
        ParamUnits={'pixels','pixels','#','#'};
        ParamDescription={'x-position', 'y-position', 'Intensity', 'background'};
        nHyperParams=5;
        HyperParamNames= {'Beta_pos', 'Mean_I', 'Kappa_I', 'Mean_bg', 'Kappa_bg'};
    end % constant properties

    properties (Access=protected)
        GPUGaussMLEFitType=1;
    end
    
    methods (Access=public)
        function obj = Gauss2DMAP(imsize_,psf_sigma_)
            % obj = Gauss2DMAP(imsize,psf_sigma) - Make a new Gauss2DMAP for
            % point localization in 2D with a fixes PSF.
            % (in) imsize: scalar int - size of image in pixels on each side (min: obj.MinSize)
            % (in) psf_sigma: scalar double>0 - size of PSF in pixels
            % (out) obj - A new object
            obj@MappelBase(@Gauss2DMAP_Iface, imsize_, psf_sigma_);
        end
 
    end %public methods
end % classdef
