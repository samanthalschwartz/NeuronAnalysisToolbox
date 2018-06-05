function  RGB=hyperCM(lambda, mask)
    % RGB = hyperCM(lambda)
    % [in] lambda - Reverse sorted vector of wavelengths to map to RGB in (nm)
    % [in] mask - [optional] 1XN logical false=keep wavelength true=remove wavelength
    % [out] rgb - Nx3 array with the RGB color values for each input lambda in each column  
    % 
    % Mark Olah (mjo@cs.unm.edu) 02-10-2014 (revised 07-15-2014)
    %
    % Values closer to min will be represented as blues and values closer to max will be 
    % represented as reds. This is intended to replace rgbstretch which does not properly span
    % the spectrum.
    %
    % Based on: http://en.wikipedia.org/wiki/File:HSV-RGB-comparison.svg
    % We leave off the 300deg-360deg part of the mapping which would bring purple
    % back to red.  Also we add a bin at long wavelengths that is constant red and a bin at 
    % short wavelengths that is constant purple.
  
    nBins=7;

    %This spaces the bins evenly where 1/nBins of lambdas are in each bin
    bins=lambda(uint64(linspace(length(lambda),1,nBins+1))); 
    bins(end)=bins(end)+1e-6; %Correct last bin endpoint to be just past the longest wavelength

    % lbins{i} is a logical the same size as lambda and is true if it is the i'th bin, 
    % which means lambda is between the i'th and i+1'st threshold barrier between bins
    lbins=cellmap(@(b,t) (b<=lambda) & (lambda<t), bins(1:end-1), bins(2:end));

    % These are the scaled lambda values for each bin giving relative distance from
    % start of bin to the end of bin
    lbvals=cellmap(@(b,t) (lambda-b)/(t-b), bins(1:end-1), bins(2:end));

    % Do the mapping of Bins to R G and B.
    R=lbins{7}+lbins{6} + lbins{5}.*lbvals{5} + lbins{2}.*(1-lbvals{2}) + lbins{1};
    G=lbins{6}.*(1-lbvals{6}) + lbins{5} + lbins{4} + lbins{3}.*lbvals{3};
    B=lbins{4}.*(1-lbvals{4}) + lbins{3} + lbins{2} + lbins{1};
    RGB=[R' G' B'];
    if nargin==2
        RGB(mask,:)=0.0;
    end
end
