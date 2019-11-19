%% ALMC Pipeline Class
classdef SRPipeline < handle
    properties
        %% define parameters
        raw_data_file = [];
        save_data_path = [];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % acquisition parameters
        acq_nchannels = 2;          % number of channels acquired in experiment
        acq_channel_shift = [0 0; 0 0; 0 0];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Temporal Median Filter parameters
        TMF_do_TMF = 'true';      % should TMF be performed? ('true' or 'false')
        TMF_medfiltradius = 51;
        TMF_keyframe_dist = 10;
        TMF_quantile = 20;
        TMF_save_TMF = 'yes';     % save the temporal median filtered data? options: 'no', 'tif', 'mat'. 'tif': saves each TMF corrected channel as a tif, 'mat': saves the TMF corrected channel in the output Matlab file
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % all the parameters that define single moclecule localization with ThunderSTORM
        Tstorm_camera_pixelsize = 100;  % camera pixelsize in nm
        Tstorm_camera_n_pixel_x =[];
        Tstorm_camera_n_pixel_y = [];
        Tstorm_camera_photons2adu = 3.78; % photoelectrons per A/D count, see camera test spec sheet
        Tstorm_camera_offset = 0;         % Base level [A/D counts], what is the lowest count that the camera is clamped to, this is 0 if TMF (temporal median filtering is used)
        Tstorm_camera_isemgain = 'true';  % is EMCCD gain being used, 'true' or 'false'
        Tstorm_camera_emgain = 100;       % EM gain value
        
        Tstorm_analysis_image_filter = 'filter=[Wavelet filter (B-Spline)] scale=2.0 order=3 ';
        % Image filter applied before any analysis
        % Options (with typical values):
        % wavelet filter:                   'filter=[Wavelet filter (B-Spline)] scale=2.0 order=3 '
        % averaging box filter:             'filter=[Averaging (Box) filter] size=3 '
        % difference of averaging filter:   'filter=[Difference of averaging filters] size1=3 size2=5 '
        % lowered gaussian filter:          'filter=[Gaussian filter] sigma=1.6 '
        % difference of gaussian filter:    'filter=[Difference-of-Gaussians filter] sigma2=1.6 sigma1=1.0 '
        % median filter:                    'filter=[Median filter] size=3 pattern=box ', or 'pattern=cross '
        % no filter:                        'filter=[No filter] '
        
        Tstorm_analysis_approx_molecule_position = 'detector=[Local maximum] connectivity=8-neighbourhood threshold=2.5*std(Wave.F1) ';
        % Determines method for initial guess of molecule localization
        % Options:
        % Local maximum:                        'detector=[Local maximum] connectivity=8-neighbourhood threshold=std(Wave.F1) ', connectivity is either '8-neighbourhood' or '4-neighbourhood'
        % Centroid of connected components:     'detector=[Centroid of connected components] watershed=false threshold=std(Wave.F1) '
        % Non-maximum suppression:              'detector=[Non-maximum suppression] threshold=std(Wave.F1) radius=1 '
        % please note: threshold can be specified. var(),std(),mean(),median(),min(),max(),sum(),abs(), normal opterators apply (+,-,/,*), F: filtered image, I: unfiltered image
        
        Tstorm_analysis_molecule_localization = 'estimator=[PSF: Integrated Gaussian] sigma=2 fitradius=5 method=[Weighted Least squares] full_image_fitting=false mfaenabled=false ';
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
        
        Tstorm_analysis_visualization = 'renderer=[No Renderer]';
        % Selects if and how the analysis results are visualized
        % Options:
        % No visusalization:                '[No Renderer]'
        % Averaged shifted histograms:      'renderer=[Averaged shifted histograms] magnification=5.0 colorizez=false threed=false shifts=2 repaint=50'
        % Scater plot:                      'renderer=[Scatter plot] magnification=5.0 colorizez=false threed=false repaint=50'
        % Normalized Gaussian:              'renderer=[Normalized Gaussian] dxforce=false magnification=5.0 dx=20.0 colorizez=false threed=false dzforce=false repaint=50'
        % Histograms:                       'renderer=Histograms magnification=5.0 avg=0 colorizez=false threed=false repaint=50'
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Drift correction parameters
        drift_onoff = 'on';                 % turn drift correction on/off, Important if drift has been corrected prior to running the pipeline, e_g. using CAT_app
        drift_method = 'RCC';               % select drift correction method, 'RCC' (redundant cross correlation), 'MCC' (mean cross correlation), 'DCC' (direct cross correlation)
        drift_segpara = 500;                % segmentation parameters, how many frames are put together
        drift_binsize = 10;                 % in nm, binsize used in cross correlation
        drift_rmax = 15;                    % error threshold for re-calculate the drift (pixel) only necessary for 'RCC'
        drift_align = 'no';                 % only necessary when two sequential channels are present.
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
        image_gen_generate_image = 'true';
        image_gen_pixelsize_b = 10;         % pixelsize of binary image in nm
        image_gen_pixelsize_g = 5;         % pixelsize of gaussian image in nm
        image_gen_precision_filter = [0 30]; % only use localizations that are within the interval [min(p.image_gen.precision_filter),max(p.image_gen.precision_filter)] (in nm)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%End of parameters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
   %% functions here ---
    methods
        %% main function to run things
        function AnalyzeDataFiles(obj)
            % this is the main method to run in order to loop through all the files in the raw_data_file
            % cell array:
            numdatafiles = size(obj.raw_data_file,1);
            %----- loop through each acquisition (multichannels)
            for ff = 1:numdatafiles 
                disp(['Experiment ', num2str(ff),' of ',num2str(numdatafiles)]);
                % loop through each ch of acquisition and build up data
                % file with fields for each channel result
                data = [];
                currfilename = [];
                concat_savefilename = [];
                %--- loop throgh each channel of the acquisition
                for cc = 1:obj.acq_nchannels 
                    disp(['Channel ',num2str(cc),' of ',num2str(obj.acq_nchannels)])
                    datapath = obj.raw_data_file{ff,(2*cc - 1)};
                    datafilename = obj.raw_data_file{ff,(2*cc)};
                    % keep filename string for building up save file name
                    % for acquisition
                    [~,save_str,~] = fileparts(datafilename);
                    if strcmp(currfilename,datafilename) % check if filenames are the same, 
                    % if so then there may be multiple channels within the image, so pass the channel info into the bf_load_parts_v2
                    % note that this will only be true if cc>1 because
                    % currfilename is not set until after first loop
                        channel = cc;
                        concat_savefilename = [concat_savefilename '_' num2str(cc)]; % only add on the channel number
                    else
                        channel = 1;
                        concat_savefilename = [concat_savefilename '_' save_str]; % add in the file name
                    end
                    currfilename = datafilename;
                    data.SML_raw_data.(['ch' num2str(cc)]) = obj.loadANDfitImage(datapath,datafilename,channel);
                end
                
                %--- drift correct if on
                if obj.drift_onoff
                    data = driftcorrect(obj,data);
                else
                    data.SML_data = data.SML_raw_data;
                end
                
                %---- make the plots here of reconstructed acquisiiton -
                %all channels
                if strcmp(obj.image_gen_generate_image,'true')
                    obj.genImage_binary(data,concat_savefilename);
                    obj.genImage_gauss(data,concat_savefilename);
                else
                    disp('No images are generated');
                end
                %---- save the results
                disp('Saving data');
                p = obj.makePstruct();
                SMLMpipeline_save_data_v2(data,p,obj.raw_data_file(ff,:),obj.save_data_path);
            end
        end  
        
        %% helper functions
        function fitresults = loadANDfitImage(obj,path,file,channel)
            % this function works on a single image file
            % loads raw data, set pixel info
            % do temporal median filtering if selected
            % run Thunderstorm and save file
            % return outputdata for further processing
            disp('Loading data');
            data=bf_load_parts_v2(path,file,channel);
            obj.Tstorm_camera_n_pixel_x = size(data,2);
            obj.Tstorm_camera_n_pixel_y = size(data,1);
            % make the individual filename for passing along and saving
            % things
            [~,NAME,~] = fileparts(file);
            chfilename = [NAME '_ch',num2str(channel)];
            % do temp med filter if it is set
            if obj.TMF_do_TMF
                [data,~] = obj.TempMedFilter(data);
                % generate TMF file save name base s .mat or .tif depending on obj.TMF_save_TMF  
                savename = obj.makeTMFsavename(chfilename); 
                if ~isempty(savename) %if savename is empty then don't save
                fullsavefilepathTMF = fullfile(path,savename);
                display(['Saving TMF image as: ' fullsavefilepathTMF]);
                obj.saveImage(data,fullsavefilepathTMF);
                end
            end
            fullsavefilenameANDpath = fullfile(path,chfilename);
            % send image to MIJI/Thunderstorm, run analysis based on class
            % params, then import back into Matlab using load_spreadsheet
            % function
            fitresults = obj.runThunderstormReturnData(data,fullsavefilenameANDpath);
            MIJ.run("Close All");
            MIJ.run("Close All");
            MIJ.run("Collect Garbage");
        end
        function get_raw_data_files(obj)
            % calls the gen_raw_data_file_variable function to generate
            % cell array of files
            % set object parameters obj.raw_data_file
            obj.raw_data_file = gen_raw_data_file_variable();
            numchannels = size(obj.raw_data_file,2)/2;
            % check that num channels is the same as what the user
            % specified in obj.acq_nchannels
            display(['Based on selected file input setting acq_nchannels to ' num2str(obj.acq_nchannels)]);
            obj.acq_nchannels = numchannels;
            
        end
        function selectSaveDir(obj)
            % prompts users to select a save directory, starts with first
            % file in list as example
            % sets object parameter obj.save_data_path
            obj.save_data_path = uigetdir(obj.raw_data_file{1,1},'Select folder where output data will be saved');
        end     
        function save(obj)
            if isempty(obj.save_data_path)
                obj.selectSaveDir();
            end
            save(fullfile(obj.save_data_path,['SRPipeline_' datestr(datetime('now'))]));
        end
        
        %% temp median filter related functions
        function  [data_tmf,background] = TempMedFilter(obj,data)
            disp('Temporal Median Filtering is running');
            [data_tmf,background]=TemporalMedianFilter3(data,obj.TMF_medfiltradius,obj.TMF_keyframe_dist,obj.TMF_quantile);
        end
        function savename = makeTMFsavename(obj,basename)
            % save the temporal median filtered data? options: 'no', 'tif', 'mat'. 'tif': saves each TMF corrected channel as a tif, 'mat': saves the TMF corrected channel in the output Matlab file
            switch obj.TMF_save_TMF
                case 'no'
                    savename = [];
                case 'yes'
                    savename = [];
                case 'mat'
                    savename = [basename '.mat'];
                case 'tif'
                    savename = [basename '.tif'];
            end
        end
        %% set camera params
        function detectANDsetCameraOffset(obj,data)
            % sets camera offset to min value of input image series
            obj.Tstorm_camera_offset = min(data(:));
        end
        %% -- functions for drift correct
        function outputdata = driftcorrect(obj,inputdata)
            disp('Drift correction');
            p = obj.makePstruct(); % reformat obj params into p structure for drift correct call
            p.acq.flag_same_file=false; % this is still incorrect to set here
            outputdata=SMLMpipeline_corr_drift_v2(p,inputdata);
        end
        function p = makePstruct(obj)
            p.acq = obj.acq_nchannels;
            p.Tstorm.camera.pixelsize = obj.Tstorm_camera_pixelsize;
            p.Tstorm.camera.photons2adu = obj.Tstorm_camera_photons2adu;
            p.Tstorm.camera.offset = obj.Tstorm_camera_offset;
            p.Tstorm.camera.isemgain =  obj.Tstorm_camera_isemgain;
            p.Tstorm.camera.emgain = obj.Tstorm_camera_emgain;       % EM gain value
            p.drift.onoff = obj.drift_onoff;
            p.drift.method = obj.drift_method;
            p.drift.segpara = obj.drift_segpara;
            p.drift.binsize = obj.drift_binsize;
            p.drift.rmax = obj.drift_rmax;
            p.drift.align = obj.drift_align;                 % only necessary when two sequential channels are present.
        end
        %% -- functions to interact with FIJI
        function runThunderstorm(obj,data)
            % load data into fiji and then run Thunderstorm using obj
            % params
            SRPipeline.loadImage2FIJI(data);
            disp('Starting ThunderSTORM analysis');
            MIJ.run('Camera setup', ['offset=',num2str(obj.Tstorm_camera_offset), ' isemgain=',obj.Tstorm_camera_isemgain, ' photons2adu=',num2str(obj.Tstorm_camera_photons2adu), ' gainem=',num2str(obj.Tstorm_camera_emgain) ,' pixelsize=', num2str(obj.Tstorm_camera_pixelsize)]);
            MIJ.run('Run analysis', [obj.Tstorm_analysis_image_filter, obj.Tstorm_analysis_approx_molecule_position, obj.Tstorm_analysis_molecule_localization, obj.Tstorm_analysis_visualization]);
        end
        
        function outputdata = runThunderstormReturnData(obj,inputdata,fullsavefilename)
            obj.runThunderstorm(inputdata);
            outputdata = SRPipeline.saveANDloadThunderstormLocalizationData(fullsavefilename);
        end
        %% -- functions to make images
        function genImage_binary(obj,data,concat_savefilename)
            disp('Generating Tiff image');
            % generate binary image
            merge_image_str = []; %string for IJ Merge Channels call: build this up as you loop through each channel
            for k=1:obj.acq_nchannels
                imagename = ['binary_ch' num2str(k)];   
                image = generate_superresimage_v2(data.SML_data.(str).position_x,...
                    data.SML_data.(str).position_y,data.SML_data.(str).precision,...
                    data.SML_data.(str).photons,ones(size(data.SML_data.(str).photons)),...
                    obj.image_gen_pixelsize_b,'binary',...
                    [obj.image_gen_pixelsize_b/2 obj.Tstorm_camera_n_pixel_x*obj.Tstorm_camera_pixelsize-obj.image_gen_pixelsize_b/2 obj.image_gen_pixelsize_b/2 obj.Tstorm_camera_n_pixel_y*obj.Tstorm_camera_pixelsize-obj.image_gen_pixelsize_b/2],...
                    'precision',obj.image_gen_precision_filter);
                MIJ.createImage(imagename,uint16(image),true);
                merge_image_str = [merge_image_str, strcat('c',num2str(k)),'=', imagename, ' '];
            end
            if obj.acq_nchannels>1
                MIJ.run("Merge Channels...", merge_image_str,' create' );
            end
            save_str = strcat(concat_savefilename,'_b',num2str(obj.image_gen_pixelsize_b),'.tif');
            MIJ.run("Properties...", ['channels=',num2str(obj.acq_nchannels),' slices=1 frames=1 unit=um pixel_width=',num2str(obj.image_gen_pixelsize_b/1000),' pixel_height=',num2str(obj.image_gen_pixelsize_b/1000),' voxel_depth=1']);
            MIJ.run("Save", ['Tiff..., path=[',fullfile(obj.save_data_path,save_str),']']);
            MIJ.run("Close All");
        end
        function genImage_gauss(obj,data,concat_savefilename)
            disp('Generating Gaussian image');
            % generate binary image
            merge_image_str = []; %string for IJ Merge Channels call: build this up as you loop through each channel
            for k=1:obj.acq_nchannels
                imagename = ['gaussian_ch' num2str(k)];   
                image = generate_superresimage_v2(data.SML_data.(str).position_x,...
                    data.SML_data.(str).position_y,data.SML_data.(str).precision,...
                    data.SML_data.(str).photons,ones(size(data.SML_data.(str).photons)),...
                    obj.image_gen_pixelsize_g,'gaussian',...
                    [obj.image_gen_pixelsize_g/2 obj.Tstorm_camera_n_pixel_x*obj.Tstorm_camera_pixelsize-obj.image_gen_pixelsize_g/2 obj.image_gen_pixelsize_g/2 obj.Tstorm_camera_n_pixel_y*obj.Tstorm_camera_pixelsize-obj.image_gen_pixelsize_g/2],...
                    'precision',obj.image_gen_precision_filter);
                MIJ.createImage(imagename,uint16(image),true);
                merge_image_str = [merge_image_str, strcat('c',num2str(k)),'=', imagename, ' '];
            end
            if obj.acq_nchannels>1
                MIJ.run("Merge Channels...", merge_image_str,' create' );
            end
            save_str = strcat(concat_savefilename,'_g',num2str(obj.image_gen_pixelsize_g),'.tif');
            MIJ.run("Properties...", ['channels=',num2str(obj.acq_nchannels),' slices=1 frames=1 unit=um pixel_width=',num2str(obj.image_gen_pixelsize_g/1000),' pixel_height=',num2str(obj.image_gen.pixelsize_g/1000),' voxel_depth=1']);
            MIJ.run("Save", ['Tiff..., path=[',fullfile(obj.save_data_path,save_str),']']);
            MIJ.run("Close All");
        end
        
        %% --- helper function to view parameters
        function viewTMFparams(obj)
            display(['TMF_do_TMF: ' obj.TMF_do_TMF]);    % should TMF be performed? ('true' or 'false')
            display(['TMF_medfiltradius: ' obj.TMF_medfiltradius]);
            display(['TMF_keyframe_dist: ' num2str(obj.TMF_keyframe_dist) ]);
            display(['TMF_quantile: ' num2str(obj.TMF_quantile)]);
            display([ 'TMF_save_TMF: ' obj.TMF_save_TMF   ]);   % save the temporal median filtered data? options: 'no', 'tif', 'mat'. 'tif': saves each TMF corrected channel as a tif, 'mat': saves the TMF corrected channel in the output Matlab file
        end
        
        function viewCameraparams(obj)
            display(['Tstorm_camera_pixelsize: ' num2str(obj.Tstorm_camera_pixelsize)]);
            display(['Tstorm_camera_n_pixel_x: ' num2str(obj.Tstorm_camera_n_pixel_x)]);
            display(['Tstorm_camera_n_pixel_y: ' num2str(obj.Tstorm_camera_n_pixel_y)]);
            display(['Tstorm_camera_photons2adu: ' num2str(obj.Tstorm_camera_photons2adu)]);
            display(['Tstorm_camera_offset: ' num2str(obj.Tstorm_camera_offset)]);
            display(['Tstorm_camera_isemgain: ' obj.Tstorm_camera_isemgain]);
            display(['Tstorm_camera_emgain: ' num2str(obj.Tstorm_camera_emgain)]);
        end
        
        function viewAnalysisparams(obj)
            display(['Tstorm_analysis_image_filter: ' obj.Tstorm_analysis_image_filter]);
            display(['Tstorm_analysis_approx_molecule_position: ' obj.Tstorm_analysis_approx_molecule_position]);
            display(['Tstorm_analysis_molecule_localization: ' obj.Tstorm_analysis_molecule_localization]);
            display(['Tstorm_analysis_visualization: ' obj.Tstorm_analysis_visualization]);
        end
        function viewdriftparams(obj)
            display(['drift_onoff: ' obj.drift_onoff]);
            display(['drift_method: ' obj.drift_method]);
            display(['drift_segpara: ' num2str(obj.drift_segpara)]);
            display(['drift_binsize: ' num2str(obj.drift_binsize)]);
            display(['drift_rmax: ' num2str(obj.drift_rmax)]);
            display(['drift_align: ' obj.drift_align]);
        end
        function viewImageGenparams(obj)
            display(['image_gen_generate_image: ' obj.image_gen_generate_image]);
            display(['image_gen_pixelsize_b: ' num2str(obj.image_gen_pixelsize_b)]);
            display(['image_gen_pixelsize_g: ' num2str(obj.image_gen_pixelsize_g)]);
            display(['image_gen_precision_filter: ' num2str(obj.image_gen_precision_filter)]);
        end
        
    end
    methods(Static)
        % methods not dependent on the class here:
%         can call them by typing "SRPipeline.startMIJI()" in command
%         prompt or "SRPipeline.loadImage2FIJI(data)" where data is an
%         image sequence. Non static methods above are dependent on parameters
%         from the class, these methods are not
        function startMIJI()
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
        end
        function loadImage2FIJI(data)
            disp('Transferring data into FIJI');
            MIJ.createImage('SMLM raw data',data,true);
        end
        function data = saveANDloadThunderstormLocalizationData(fullsavefilenameANDpath)
            % this saves thunderstorm data open in MIJI to csv file and
            % reloads into matlab
            %filename = ['\results_file' num2str(j) '_ch_' num2str(i) '.csv'];
            disp('ThunderSTORM analysis is finished, transferring SML data into Matlab');
            MIJ.run("Export results", ['filepath=[' fullsavefilenameANDpath, '] fileformat=[CSV (comma separated)] sigma=true intensity=true chi2=true offset=true saveprotocol=true x=true y=true bkgstd=true id=true uncertainty=true frame=true']);
            [FILEPATH,NAME,EXT] = fileparts(fullsavefilenameANDpath);
            data = load_spreadsheet([FILEPATH,filesep],[NAME,EXT],'ThunderSTORM');
            %         delete(strcat(save_data_path,'dmy.csv'));
            %         delete(strcat(save_data_path,'dmy-protocol.txt'));
            MIJ.closeAllWindows %closes all open windows in the current Miji instance, necessary to not run into Java heap space issues
        end
        function saveImage(image,fullsavefilename)
            options.overwrite = true;
            saveastiff(image,fullsavefilename,options)
        end
    end
end
%---