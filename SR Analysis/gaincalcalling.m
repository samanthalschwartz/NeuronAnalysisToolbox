namestr = '250';
gain = loadtiff(['F:\SR\GainCorrect\180724_gain' namestr '.tif']);
% gain = cat(3,gain1{1},gain1{2});

bg = loadtiff(['F:\SR\GainCorrect\180724_bg' namestr '.tif']);
% fullbg = cat(3,bg{1},bg{2});

out = cal_readnoise(gain,bg,100,-1,0);
saveas(gcf,['F:\SR\GainCorrect\gain' namestr 'cresults'],'png');
