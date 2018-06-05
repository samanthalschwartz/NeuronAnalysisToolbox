function testEmanuel()
%     emanuel_root_dir = '/olah/home/mjo/LidkeLab/Data/Emanuel/';
    emanuel_root_dir = '/home/mjo/LidkeLab/Data/Emanuel/';
%     emanuel_root_dir = 'C:\Users\mjo\LidkeLab\Data\Emanuel';
    emanuel_data_dir = fullfile(emanuel_root_dir,  '150624_A431_EGF-QD585_EGF-QD655');
    
    cr1=fullfile( emanuel_data_dir, 'ChannelRegistation-2015-6-24-14-7-6.mat');
    cr2=fullfile( emanuel_data_dir, 'ChannelRegistation-2015-6-24-15-48-26.mat');
    rmse = RegisteredPairAnalysis.testChannelRegistrationV1(cr1,cr2);
    
    rpt_dir = fullfile(emanuel_data_dir,'RPT');
    rpt1 = fullfile(rpt_dir, 'cell8_welll3_160pM-EGF-QD-each_set4-2#0001-2015-6-24-15-22-15_LEFT1.rpt');
    rpt2 = fullfile(rpt_dir, 'cell8_welll3_160pM-EGF-QD-each_set4-2#0001-2015-6-24-15-22-15_RIGHT1.rpt');
    disp('Tform RMSE');
    disp(rmse);
    
    pa = RegisteredPairAnalysis(rpt1, rpt2, cr1, [585, 655]);
    pa.plotTracks();
    icp = pa.makeInteractionChangePoint( [1,2]);
    icp.plotPairMatrixDistances();
%     bx = Boxxer2D(
    
end
