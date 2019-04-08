function image_out = shift_withpadding(image_in,shiftvec,padwidth)
padded_image_in = padzeros(image_in,padwidth);
shiftim = shift(padded_image_in,shiftvec,1);
image_out = shiftim(padwidth:(end-padwidth),padwidth:(end-padwidth));
end