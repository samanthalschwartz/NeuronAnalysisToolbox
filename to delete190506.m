test = reshape(image,[size(image,1),size(image,2),15,27]);
testm = squeeze(max(test,[],3));
dipshow(testm,'percentile')
uiopen('Y:\Sam\190503 GabaA2-HaloTag\FRAP\Cry2Olig_HaloTag647_20190506_91406 AM\Cry2Olig_HaloTag647_t0002_w0001.tif',1)