function img_out = LoG_threshold(im_in,gsig,lsig)
gim = gaussf(im_in,gsig);
lim = dxx(gim,lsig) + dyy(gim,lsig);
lim = -lim;
lim(lim<0) = 0;
img_out = lim;
end