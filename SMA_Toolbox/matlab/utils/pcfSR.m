%pcfSR Calculated the photon conversion factor of poisson noise
%    from the image data only
%
% SYNOPSIS:
%   out = pcf(in, kcut, gN, plotFlag)
%
%   kcut:frequencies that should be cut away
%        in the high pass filtering [0, sqrt(dimension of in)]
%        the 'corners' of the rectangular image have freq > 1
%
% DEFAULTS:
%   kcut = 1, gN = 0, plotFlag = 0;
%
% SEE ALSO
%   noise with 'poisson' option and paramter
%
% EXAMPLE:
%   a = readim;
%   psf = gaussf(deltaim,2);%not really a point spread function of a microscope
%   ap = convolve(a,psf);
%   apn = noise(ap,'poisson',3.8,0);
%   pcf(apn)
%
% LITERATURE:
%   R. Heintzmann, B. Rieger, G. Ficz, T.M. Jovin, Co-localization of noisy data, submitted

% (C) Copyright 2004           Department of Molecular Biology
%     All rights reserved      Max-Planck-Institute for Biophysical Chemistry
%                              Am Fassberg 11, 37077 G"ottingen
%                              Germany
%
% Bernd Rieger, Rainer Heintzmann August 2004.
% March 2015, and read noise and remove axis in FT as by Peter Relich (BR)
% April 2015, modified function for Super Resolution Time Stacked Data (PR)


function [gain, offset] = pcfSR(varargin)
d = struct('menu','060',...
    'display','Estimate photon conversion factor',...
    'inparams',struct('name',       {'in','kcut','gN','plotFlag'},...
    'description',{'Input image','Cut frequencies larger than','Read noise','Plot Flag'},...
    'type',       {'image','array','array','array'},...
    'dim_check',  {0,-1,-1,-1},...
    'range_check',{[],'R+','R+','R+'},...
    'required',   {1,0,0,0},...
    'default',    {'a',1,0,0}...
    ),...
    'outparams',struct('name',{'gain','offset'},...
    'description',{'Gain','Offset'},...
    'type',{'image','image'}...
    )...
    );
try
    [in,kcut,gN,plotFlag] = getparams(d,varargin{:});
catch
    if ~isempty(paramerror)
        error(paramerror)
    else
        error(firsterr)
    end
end

if kcut <0 || kcut >sqrt(ndims(in))
    error('kcut in [0,sqrt(dimensions of in)]');
end
% modified for 3D (PKR)
if ndims(in)~=3
    error('assuming 3D (x,y,t) input');
end

% explicitly modifying this step to work with time stacks (PKR)
sz = imsize(in);
if sz(1)~=sz(2)
    fprintf('Cutting the image to a square size of %d\n',min(sz(1),sz(2)));
    in = cut(in,[min(sz),min(sz),sz(3)]);
    sz = imsize(in);
end

% parse out frames in the stack to analyze (hard limit to 100 for now)
if sz(3)>101 % 100 somewhat evenly spaced frames if there are more than 100
   ratio = (sz(3)-1)/100;
   vecind = ceil((1:100)*ratio);
end
if sz(3)<= 100 % can only parse as many frames are available
   vecind = 1:sz(3)-1; % never pick frame 0
end

% pre-initialize the fourier transform image stack and variance vector
f = newim(sz(1),sz(2),length(vecind),'scomplex');
Neff = zeros(length(vecind),1);
TotalVar = zeros(length(vecind),1);
SumPix = zeros(length(vecind),1);

% constant values through all images in the stack
m = rr(sz(1),sz(2))>(kcut*size(in,1)/2);
ssz = imsize(m);
sz2 = ssz./2;
% spectral cross cut off
n = 1;
m(floor(sz2(1))-n:ceil(sz2(1))+n,:)=0;
m(:,floor(sz2(2))-n:ceil(sz2(2))+n)=0;
volf = mean(m);%sum(m)/prod(size(m))

% calculate the variance at all sampled points
for ii = 1:length(vecind)
    f(:,:,ii-1) = ft(squeeze(in(:,:,vecind(ii))));
    nf = m*squeeze(f(:,:,ii-1));
    Neff(ii) = sum(abs(nf)^2)-gN^2*sum(m);
    TotalVar(ii) = Neff(ii)/volf;
    SumPix(ii) = sum(in(:,:,vecind(ii)));
end

% loop over the image slices to generate the statistics
[tempoffset, stats] = robustfit(SumPix, TotalVar);
offset = -tempoffset(1)/(ssz(1)^2)/ tempoffset(2);
%SDoffset = stats.se(1)/(horsize*vertsize)/gain;
Allvar = sum(TotalVar);
gain = Allvar/sum(in(:,:,vecind)-offset);

if plotFlag
    figure;
    plot(SumPix,TotalVar,'+')
    hold on
    title('Gain Offset Regression');
    xlabel('Sum Image (ADUs)');
    ylabel('Total Variance (ADU^2s)');
    set(gca,'FontSize',12,'FontWeight','bold');
    plot(min(SumPix):(max(SumPix)-min(SumPix))/100:max(SumPix),...
        tempoffset(1)+tempoffset(2)*(min(SumPix):(max(SumPix)-min(SumPix))/100:max(SumPix)),'r')
    hold off
end



