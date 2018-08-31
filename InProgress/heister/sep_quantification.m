uiopen('E:\Matt Becker Data (For Review)\SEPGlua1_mch\20150609_SEPGlua1_mch_4_2.tif',1)
ga.imgDggCutoff(gaussf(sep))
image_out = varif(ans,7,'elliptic')
ga.imgThreshold_fixedUserInput(image_out)
spines = ans;
ga.viewMaskOverlayPerimStatic(sep,spines)