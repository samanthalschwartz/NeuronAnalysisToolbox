function sim=sumImage2D(im)
    % This is what makes the sumImage representations in the SPData and related classes.
    % It is a combination of the max and the mean images that seems to give acceptable results for most
    % data types.
    % [IN]
    % im - 3D image size: [X Y Z] any numeric type
    % [OUT]
    % sim - 2D sum/mean combination image size: [X Y] type: single

    im=single(im); %ensure image is treated like a single and resulting sim will also be single
    max_im=max(cosmicNorm(im),[],3);
    mean_im=single(mean(im,3));
    mean_im=mean_im./max(mean_im(:));
    blend=0.8;
    sim=blend*max_im+(1-blend)*mean_im;
end
