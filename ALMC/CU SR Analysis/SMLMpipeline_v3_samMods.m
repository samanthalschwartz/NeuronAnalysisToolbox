% pipeline for SMLM data analysis
% in order to run this routine you need FIJI, the ThunderSTORM plugin and
% MIJ. Follow the instructions on http://bigwww.epfl.ch/sage/soft/mij/
% You do not have to call MIJ or MIJ.start before calling this file.

%% define parameters
p = struct;     %contains all the parameters
data = struct;  %contains all the data 
% acquisition parameters
p.acq.nchannels = 2;          % number of channels acquired in experiment
p.acq.channel_shift = [0 0; 0 0; 0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temporal Median Filter parameters 
p.TMF.do_TMF = 'true';      % should TMF be performed? ('true' or 'false')
p.TMF.medfiltradius = 51;
p.TMF.keyframe_dist = 10;
p.TMF.quantile = 20;
p.TMF.save_TMF = 'yes';     % save the temporal median filtered data? options: 'no', 'tif', 'mat'. 'tif': saves each TMF corrected channel as a tif, 'mat': saves the TMF corrected channel in the output Matlab file 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% all the parameters that define single moclecule localization with ThunderSTORM
p.Tstorm.camera.pixelsize = 100;  % camera pixelsize in nm
p.Tstorm.camera.photons2adu = 3.78; % photoelectrons per A/D count, see camera test spec sheet
p.Tstorm.camera.offset = 0;         % Base level [A/D counts], what is the lowest count that the camera is clamped to, this is 0 if TMF (temporal median filtering is used)
p.Tstorm.camera.isemgain = 'true';  % is EMCCD gain being used, 'true' or 'false'
p.Tstorm.camera.emgain = 100;       % EM gain value

p.Tstorm.analysis.image_filter = 'filter=[Wavelet filter (B-Spline)] scale=2.0 order=3 '; 
% Image filter applied before any analysis 
% Options (with typical values):
% wavelet filter:                   'filter=[Wavelet filter (B-Spline)] scale=2.0 order=3 '
% averaging box filter:             'filter=[Averaging (Box) filter] size=3 '
% difference of averaging filter:   'filter=[Difference of averaging filters] size1=3 size2=5 '
% lowered gaussian filter:          'filter=[Gaussian filter] sigma=1.6 '
% difference of gaussian filter:    'filter=[Difference-of-Gaussians filter] sigma2=1.6 sigma1=1.0 '
% median filter:                    'filter=[Median filter] size=3 pattern=box ', or 'pattern=cross '
% no filter:                        'filter=[No filter] '

p.Tstorm.analysis.approx_molecule_position = 'detector=[Local maximum] connectivity=8-neighbourhood threshold=2.5*std(Wave.F1) ';
% Determines method for initial guess of molecule localization
% Options:
% Local maximum:                        'detector=[Local maximum] connectivity=8-neighbourhood threshold=std(Wave.F1) ', connectivity is either '8-neighbourhood' or '4-neighbourhood'
% Centroid of connected components:     'detector=[Centroid of connected components] watershed=false threshold=std(Wave.F1) '
% Non-maximum suppression:              'detector=[Non-maximum suppression] threshold=std(Wave.F1) radius=1 '
% please note: threshold can be specified. var(),std(),mean(),median(),min(),max(),sum(),abs(), normal opterators apply (+,-,/,*), F: filtered image, I: unfiltered image

p.Tstorm.analysis.molecule_localization = 'estimator=[PSF: Integrated Gaussian] sigma=2 fitradius=5 method=[Weighted Least squares] full_image_fitting=false mfaenabled=false ';
% Select method for singele molecule localization
% Options:
% PSF: Integrated Gaussian:                 'estimator=[PSF: Integrated Gaussian] sigma=1.5 fitradius=3 method=[Weighted Least squares] full_image_fitting=false mfaenabled=false '
% methods: [Weighted Least squares],[Least squares],[Maximum likelihood]
% if multi-emitter fit analysis is enabled then more parameters are needed: 'mfaenabled=true keep_same_intensity=false nmax=5 fixed_intensity=true expected_intensity=500:2500 pvalue=1.0E-6 '
% PSF: Gaussian:                            'estimator=[PSF: Gaussian] sigma=1.5 fitradius=3 method=[Weighted Least squares] full_image_fitting=false mfaenabled=false '    
% PSF: Elliptical Gaussian (3D astigmatism):
% Radial symmetry:                          'estimator=[Radial symmetry] fitradius=3 '
% Centroid of local neighborhood:           'estimator=[Centroid of local neighborhood] fitradius=3 '
% No estimator:                             'estimator=[No estimator] '

p.Tstorm.analysis.visualization = 'renderer=[No Renderer]';
% Selects if and how the analysis results are visualized
% Options:
% No visusalization:                '[No Renderer]'
% Averaged shifted histograms:      'renderer=[Averaged shifted histograms] magnification=5.0 colorizez=false threed=false shifts=2 repaint=50'
% Scater plot:                      'renderer=[Scatter plot] magnification=5.0 colorizez=false threed=false repaint=50'
% Normalized Gaussian:              'renderer=[Normalized Gaussian] dxforce=false magnification=5.0 dx=20.0 colorizez=false threed=false dzforce=false repaint=50'
% Histograms:                       'renderer=Histograms magnification=5.0 avg=0 colorizez=false threed=false repaint=50'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Drift correction parameters
p.drift.onoff = 'on';                 % turn drift correction on/off, Important if drift has been corrected prior to running the pipeline, e.g. using CAT_app
p.drift.method = 'RCC';               % select drift correction method, 'RCC' (redundant cross correlation), 'MCC' (mean cross correlation), 'DCC' (direct cross correlation)
p.drift.segpara = 500;                % segmentation parameters, how many frames are put together
p.drift.binsize = 10;                 % in nm, binsize used in cross correlation
p.drift.rmax = 15;                    % error threshold for re-calculate the drift (pixel) only necessary for 'RCC'
p.drift.align = 'no';                 % only necessary when two sequential channels are present. 
                                      % 'yes' will try to calculate the
                                      % drift between ch1 and ch2 that
                                      % occured between the acquisition
                                      % from background images, which is
                                      % quite unreliable. 
                                      % 'no' : suggested option! In this
                                      % case it is assumed that the
                                      % sample didn't drift between the
                                      % stop of acquistion 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters for automatically generated images
p.image_gen.generate_image = 'true';
p.image_gen.pixelsize_b = 10;         % pixelsize of binary image in nm
p.image_gen.pixelsize_g = 5;         % pixelsize of gaussian image in nm 
p.image_gen.precision_filter = [0 30]; % only use localizations that are within the interval [min(p.image_gen.precision_filter),max(p.image_gen.precision_filter)] (in nm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%End of parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% specify raw data files
if exist('raw_data_file')==1 %check first if it already exists if not start file dialog to specify data file(s)
    if iscell(raw_data_file)
        
    else
        raw_data_file={};                    % generate cell structure that contains the information which data should be analyzed
        [raw_data_file{1,2},raw_data_file{1,1}] = uigetfile('*.*', 'Load Channel 1 (recorded first)');
        
        if p.acq.nchannels == 2
            [raw_data_file{1,4},raw_data_file{1,3}] = uigetfile('*.*', 'Load Channel 2 (recorded second)',raw_data_file{1,1});
        elseif p.acq.nchannels == 3
            [raw_data_file{1,4},raw_data_file{1,3}] = uigetfile('*.*', 'Load Channel 2 (recorded second)',raw_data_file{1,1});
            [raw_data_file{1,6},raw_data_file{1,5}] = uigetfile('*.*', 'Load Channel 3 (recorded second)',raw_data_file{1,1});
        end
    end
else
        raw_data_file={};                    % generate cell structure that contains the information which data should be analyzed
        [raw_data_file{1,2},raw_data_file{1,1}] = uigetfile('*.*', 'Load Channel 1 (recorded first)');
        
        if p.acq.nchannels == 2
            [raw_data_file{1,4},raw_data_file{1,3}] = uigetfile('*.*', 'Load Channel 2 (recorded second)',raw_data_file{1,1});
        elseif p.acq.nchannels == 3
            [raw_data_file{1,4},raw_data_file{1,3}] = uigetfile('*.*', 'Load Channel 2 (recorded second)',raw_data_file{1,1});
            [raw_data_file{1,6},raw_data_file{1,5}] = uigetfile('*.*', 'Load Channel 3 (recorded second)',raw_data_file{1,1});
        end
end
%%  specify folder where data will be saved
save_data_path = uigetdir(raw_data_file{1,1},'Select folder where output data will be saved');

    %% set up MIJI
    disp('Starting FIJI');
    currentFolder=pwd;
    addpath ('C:\FIJI\scripts');
    javaaddpath 'C:\Program Files\MATLAB\R2018a\java\mij.jar';           % extend the java path to the mij.jar file
    javaaddpath 'C:\Program Files\MATLAB\R2018a\java\ij-1.52h.jar';      % extend the java path to the current ij-x.xxl.jar file
    
    %% also, make sure that the scripts/ directory of your Fiji.app/ is in MATLAB's search patch, via File  ? Set Path... (on Mac, the file chooser doesn't let you choose directories within .app packages, so you have to use the MATLAB command addpath('/Applications/Fiji.app/scripts') ). 
    Miji;
    disp('Running FIJI garbage collect');
    MIJ.run("Collect Garbage");
    cd(currentFolder); %the previous commands change the current folder, so change it back to what we started with

%% for loop that loops through all the specified experiments
%  for this to work the "specify raw data files section has to be commented
%  out and the raw data files are specified in a cell structure called 
%  'raw_data_file' where raw_data_file{i,1}= 'path', raw_data_file{i,2}= 'file1'
%  if 2 imaging channels are present then raw_data_file{i,3}= 'path' and 
%  raw_data_file{i,4}= 'file2'. If file1 == file2 then it assumes an interleaved
%  acquisition with 2 channels
for j=1: size(raw_data_file,1)
    %% for loop that loops through the channels of an acquisition
    for i=1:p.acq.nchannels
            disp(['Channel ',num2str(i),' of ',num2str(p.acq.nchannels),'   Experiment ', num2str(j),' of ',num2str(size(raw_data_file,1))]);
            disp('Loading data');
        %% load the raw data
        if i==1            
            data.raw_data=bf_load_parts_v2(raw_data_file{j,1},raw_data_file{j,2},1);
            p.Tstorm.camera.n_pixel_x = size(data.raw_data,2);
            p.Tstorm.camera.n_pixel_y = size(data.raw_data,1);
        else
            p.acq.flag_same_file = true;
            for cc = 2: p.acq.nchannels
                p.acq.flag_same_file = strcmp(raw_data_file{j,2},raw_data_file{j,2*cc}) && p.acq.flag_same_file;
            end

            if p.acq.flag_same_file==true
                data.raw_data=bf_load_parts_v2(raw_data_file{j,1},raw_data_file{j,2*i},i);
            else
                data.raw_data=bf_load_parts_v2(raw_data_file{j,2*i-1},raw_data_file{j,2*i},1);
            end
        end
        
        %% run Temporal Median Filter on the data
        if strcmp(p.TMF.do_TMF,'true')
            disp('Temporal Median Filtering is running');
            [data.data_tmf,data.background]=TemporalMedianFilter3(data.raw_data,p.TMF.medfiltradius,p.TMF.keyframe_dist,p.TMF.quantile);
%             if p.acq.nchannels==2
%                 % generate background data that can be used for aligning the
%                 % two sequantially imaged channels
%                 if ~strcmp(raw_data_file{1,2},raw_data_file{1,4})
%                     %if strcmp(p.drift.align,'yes') 
%                         data.bg_align_ch.indizes = 1 : ceil(size(data.data_tmf,3)/100) : size(data.data_tmf,3);
%                         data.bg_align_ch.indizes(end) = size(data.data_tmf,3);
%                         data.bg_align_ch.(['ch' num2str(i)]) = data.background(:,:,data.bg_align_ch.indizes);  
%                     %end                  
%                 end
%             end
%             data = rmfield(data,'background');
        else
            disp('Temporal Median Filtering is off');
        end
        
        %% load image stack into FIJI
        disp('Transferring data into FIJI');
        if strcmp(p.TMF.do_TMF,'true')
            p.Tstorm.camera.offset = min(min(min(data.data_tmf)));  % detect camera offset
            MIJ.createImage('SMLM raw data', data.data_tmf, true);
        else
            p.Tstorm.camera.offset = min(min(min(data.raw_data)));  % detect camera offset
            MIJ.createImage('SMLM raw data', data.raw_data, true);
        end

        %% run ThunderSTORM analysis
        disp('Starting ThunderSTORM analysis');
        MIJ.run('Camera setup', ['offset=',num2str(p.Tstorm.camera.offset), ' isemgain=',p.Tstorm.camera.isemgain, ' photons2adu=',num2str(p.Tstorm.camera.photons2adu), ' gainem=',num2str(p.Tstorm.camera.emgain) ,' pixelsize=', num2str(p.Tstorm.camera.pixelsize)]);
        MIJ.run('Run analysis', [p.Tstorm.analysis.image_filter, p.Tstorm.analysis.approx_molecule_position, p.Tstorm.analysis.molecule_localization, p.Tstorm.analysis.visualization]);

        %% get ThunderSTORM localization data (save as csv file first and then load into matlab from there since I haven't figured out yet how to get it directly)
        filename = ['\results_file' num2str(j) '_ch_' num2str(i) '.csv'];
        disp('ThunderSTORM analysis is finished, transferring SML data into Matlab');
        MIJ.run("Export results", ['filepath=[' fullfile(save_data_path,filename), '] fileformat=[CSV (comma separated)] sigma=true intensity=true chi2=true offset=true saveprotocol=true x=true y=true bkgstd=true id=true uncertainty=true frame=true']);
        data.SML_raw_data.(['ch' num2str(i)])=load_spreadsheet(save_data_path,filename,'ThunderSTORM');
%         delete(strcat(save_data_path,'dmy.csv'));
%         delete(strcat(save_data_path,'dmy-protocol.txt'));
        MIJ.closeAllWindows %closes all open windows in the current Miji instance, necessary to not run into Java heap space issues

        %% clean up
        data = rmfield(data,'raw_data'); 
        if strcmp(p.TMF.do_TMF,'true')
            if strcmp(p.TMF.save_TMF,'tif')
                disp('Saving temporal median filtered data to tif.');
                [~,f,~]=fileparts(raw_data_file{j,2*i});
                options.overwrite = true; 
                saveastiff(data.data_tmf, strcat(save_data_path,filesep,f,'_ch',num2str(i),'.tif'), options);
            elseif strcmp(p.TMF.save_TMF,'mat')
                data.data_TMF.(['ch' num2str(i)])=data.data_tmf;    % save temporal median filtered data in Matlab file
            end
        end
        if strcmp(p.TMF.do_TMF,'true');  data = rmfield(data,'data_tmf'); end
        MIJ.run("Close All");
        MIJ.run("Collect Garbage");
    end


    %% run drift correction algorithm
    if strcmp(p.drift.onoff,'on')   
        disp('Drift correction');
        data=SMLMpipeline_corr_drift_v2(p,data);
    else
        disp('Drift correction is turned off')
        data.SML_data = data.SML_raw_data;
    end
    %% Visualize results
    if strcmp(p.image_gen.generate_image,'true')
        disp('Generating Tiff image');
        % generate binary image
        merge_image_str = '';
        [~,save_str,~] = fileparts(raw_data_file{j,2});
        for k=1:p.acq.nchannels
            str = ['ch',num2str(k)];
            image = generate_superresimage_v2(data.SML_data.(str).position_x,data.SML_data.(str).position_y,data.SML_data.(str).precision,data.SML_data.(str).photons,ones(size(data.SML_data.(str).photons)),p.image_gen.pixelsize_b,'binary',[p.image_gen.pixelsize_b/2 p.Tstorm.camera.n_pixel_x*p.Tstorm.camera.pixelsize-p.image_gen.pixelsize_b/2 p.image_gen.pixelsize_b/2 p.Tstorm.camera.n_pixel_y*p.Tstorm.camera.pixelsize-p.image_gen.pixelsize_b/2],'precision',p.image_gen.precision_filter);
            MIJ.createImage(strcat('binary_',str),uint16(image),true);
            merge_image_str = [merge_image_str, strcat('c',num2str(k)),'=binary_',str,' '];
            if p.acq.flag_same_file == true
                save_str=strcat(save_str,str,'_');
            else
                if k>1
                    [~,dmy,~] = fileparts(raw_data_file{j,k*2}); 
                    save_str = strcat(save_str,dmy,'_');
                end
            end         
        end
        if p.acq.nchannels>1
            merge_image_str = strcat(merge_image_str,' create');
           MIJ.run("Merge Channels...", merge_image_str);
        end
        save_str = strcat(save_str,'b',num2str(p.image_gen.pixelsize_b),'.tif');
        MIJ.run("Properties...", ['channels=',num2str(p.acq.nchannels),' slices=1 frames=1 unit=um pixel_width=',num2str(p.image_gen.pixelsize_b/1000),' pixel_height=',num2str(p.image_gen.pixelsize_b/1000),' voxel_depth=1']);
        MIJ.run("Save", ['Tiff..., path=[',save_data_path,'\',save_str,']']);
        MIJ.run("Close All");
        
        % generate gaussian image
        merge_image_str = '';
        [~,save_str,~] = fileparts(raw_data_file{j,2});
        for k=1:p.acq.nchannels
            str = ['ch',num2str(k)];
           %image = generate_superresimage_v2(data.SML_data.(str).position_x,data.SML_data.(str).position_y,data.SML_data.(str).precision,data.SML_data.(str).photons,ones(size(data.SML_data.(str).photons)),p.image_gen.pixelsize_b,'binary',[p.image_gen.pixelsize_b/2 p.Tstorm.camera.n_pixel_x*p.Tstorm.camera.pixelsize-p.image_gen.pixelsize_b/2 p.image_gen.pixelsize_b/2 p.Tstorm.camera.n_pixel_y*p.Tstorm.camera.pixelsize-p.image_gen.pixelsize_b/2],'precision',p.image_gen.precision_filter);
            image = generate_superresimage_v2(data.SML_data.(str).position_x,data.SML_data.(str).position_y,data.SML_data.(str).precision,data.SML_data.(str).photons,ones(size(data.SML_data.(str).photons)),p.image_gen.pixelsize_g,'gaussian',[p.image_gen.pixelsize_g/2 p.Tstorm.camera.n_pixel_x*p.Tstorm.camera.pixelsize-p.image_gen.pixelsize_g/2 p.image_gen.pixelsize_g/2 p.Tstorm.camera.n_pixel_y*p.Tstorm.camera.pixelsize-p.image_gen.pixelsize_g/2],'precision',p.image_gen.precision_filter);           
            MIJ.createImage(strcat('gaussian_',str),uint16(image),true);
            merge_image_str = [merge_image_str, strcat('c',num2str(k)),'=gaussian_',str,' '];
            if p.acq.flag_same_file == true
                save_str=strcat(save_str,str,'_');
            else
                if k>1
                    [~,dmy,~] = fileparts(raw_data_file{j,k*2}); 
                    save_str = strcat(save_str,dmy,'_');
                end
            end 
        end
        if p.acq.nchannels>1
            merge_image_str = strcat(merge_image_str,' create');       
            MIJ.run("Merge Channels...", merge_image_str);
        end
        save_str = strcat(save_str,'g',num2str(p.image_gen.pixelsize_g),'.tif');
        MIJ.run("Properties...", ['channels=',num2str(p.acq.nchannels),' slices=1 frames=1 unit=um pixel_width=',num2str(p.image_gen.pixelsize_g/1000),' pixel_height=',num2str(p.image_gen.pixelsize_g/1000),' voxel_depth=1']);
        MIJ.run("Save", ['Tiff..., path=[',save_data_path,'\',save_str,']']);
        MIJ.run("Close All");
        
    else
        disp('No images are generated');
    end
    %% save the data
    disp('Saving data');
    SMLMpipeline_save_data_v2(data,p,raw_data_file(j,:),save_data_path);
    
end

%% close FIJI/MIJI
MIJ.exit;
cd(currentFolder);

%% clear variables from memory
clearvars data i j p save_data_path currentFolder ans options f cc dmy image k merge_image_str save_str str %raw_data_file