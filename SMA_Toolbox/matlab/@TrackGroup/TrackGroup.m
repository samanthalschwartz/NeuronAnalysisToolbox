classdef TrackGroup < Pickle

% TrackGroup - Represents a group of tracks which were identified in the
% same experimental region and can be expected to potentially interact with each other.
%
% TrackGroup is initialized with a SPTHSI object (or registered 2-channel RPT),
% and stores track localizations in two data structures.  The tracks are stored in obj.T
% which is a cell array with
% a matrix for each track.  The track matrix is described in detail below, but in general
% contains a row for each localization that comprises the track.  
%
% Statistics about each track can be gathered with methods that begin with "track", e.g.,
% obj.trackLength()
% obj.trackDuration()
% obj.trackTimeSpans()
% obj.trackSpectraBounds()
%
% Once initialized an TrackGroup can be filtered, removing tracks that do not meet
% certain quality constraints.  Track filtering is done with methods that begin with
% "filter", e.g.,
%
% filterLength
% filterDuration
% filterSpectralDeviation
% filterSelection
%
%
%
% Details:
%
% * Tracks are stored in order of length (longest to shortest)
% * Every localization is in exactly one track. And every localization from
%   every track is in the L matrix.  Thus the number of rows in L is equal
%   to the sum of all rows in all matrices of T.
% * Units: We use common units thoughout and convert to the same units from any input type
%     * Time (s)
%     * Displacments (microns)
%     * wavelength (nm)
%     * Errors are always in terms of variance (sigma^2) not standard deviations (sigma)
%
% Design decisions:
% * Really T should be a cell array of tables, but things are just too slow with
%   tables right now.  Building a cell array of tables for 4600 tracks takes
%   13.75s vs 0.27s for a cell array for matricies.  Building Localizations out
%   of T takes 3s with a table and 0.003s with an array.  Thus unfortunately we
%   will stick with matricies rather than tables internally, but will allow the
%   export of tracks as tables, cause they really are much nicer.  We store the column
%   information that would be included in the table instead as constant member variables
%   (TCol, TUnit, TDesc)
%
%
% Author: Mark J. Olah (mjo@cs.unm.edu)
% Date: Feb 2014, July 2015
%

properties (Constant = true)
        %Inherited from Pickle
        saveFileExt='.tg'; % The extension used for saved files
        SaveableDataFormats={'*.tg', 'TrackGroup (.tg)'};
        LoadableDataFormats={'*.tg;*.hsrpt;*.tsa;','All Loadable Sources (.tg, .hsrpt, .tsa)';...
                             '*.tg','TrackGroup file (.tg)';...
                             '*.tsa','TrackSegmentAnalysis file (.tsa)';...
                             '*.hsrpt', 'HSRPT object files (.hsrpt)'};
end

properties 
    frameT;    % The duration of a single frame
    ROIPhysical; % [xmin xmax ymin ymax] - The ROI bounds in physical units
    frameBounds; % [minFrameIdx, maxFrameIdx]
    Tracks; %Stored in RPT format as a cell array of matricies (uses physical units).  See RPT for details

    certainty=0.05; %The certainty to which comparisons and calculations are performed.
    
    Paths=struct(... %All filenames are relative paths from the workingDir
             'saveFile',[],... %filename(w/extension) of the .tg file where this will save to
             'data',[]... %relative path to the .spdata file we are associated with
             );
end 

properties (Transient=true)
    data;        % The SPData or HSData file associated with this track group
end


properties (Dependent=true)
    %nT: scalar - number of tracks (T)
    nT;
    %L: A matrix of all localizations for all tracks in this group
    Localizations;
    %nL: The count of all localizations for all tracks in this group
    nLocalizations;
end %dependent public properties

properties (Hidden=true)
    version=1; %For future file format version changes
end    

methods
    function obj=TrackGroup( varargin )
        % This constructor calls the obj.load method with the arguments
        % given.  Empty input arguments give an empty object, otherwise the load
        % method determines the valid input arguments
        if nargin>0
            obj.load(varargin{:});
        end
    end
    
    function load(obj, varargin)
        %
        % [Case1]: Load from saved .tg file
        %    obj.load(tsaFilePath)
        %  tsaFilePath - Full path to a .tsa file to load
        %
        % [Case2]: New TSA from SPData and RPT/SPT file
        %    obj.load(spd,roi,spt)
        %
        %  spd - An SPData filepath or object
        %  roi - An roi index into the SPData object's list of ROI,
        %        set to empty for whole area.
        %  spt/rpt - [optional] An SPT/RPT filepath or object with tracks already computed.
        %        If omitted automatically look for .rpt then .spt file when use the first save .spt file in the ROI path.
        %
        % [Case3]: New TSA from RPT
        %  obj.load(rpt)
        %  
        %  rpt - An RPT object or full path to file
        % 
        in=varargin{1};
        if ischar(in)
            [~,~,ext]=fileparts(in);
            switch ext
                case '.tsa'
                    obj.loadTSA(in);                        
                case '.spdata'
                    obj.loadSPData(varargin{:});
                case '.rpt'
                    obj.loadRPT(in);
                otherwise
                    error('TrackSegmentAnalysis:load','Unknown file format extenstion "%s"',ext);
            end
        elseif isa(in,'SPData')
            obj.loadSPData(varargin{:});
        elseif isa(in,'RPT')
            obj.loadRPT(in);
        else
            error('TrackSegmentAnalysis:load','Unknown argument type "%s"',class(in));
        end
    end
        
    function obj=TrackGroup(spt_obj, certainty)
        % TrackGroup(spt_obj) - Create a new track group from spt_object
        % (in) spt_obj - Must be a SPTHSI for now
        % TODO: Allow registered 2-channel SPT data
        if nargin>=1
            if isa(spt_obj,'SPTHSI')
                obj.initialize_SPTHSI(spt_obj);
            end
        end
        if nargin==2
            obj.set_certainty(certainty);
        else
            obj.set_certainty(0.95);
        end
    end

    function nT=get.nT(obj)
        nT=size(obj.T,1);
    end
    
    function L=get.Localizations(obj)
        L=cell2mat(obj.T);
    end

    function nL=get.nLocalizations(obj)
        nL=sum(obj.trackLengths());
    end

    function set_certainty(obj,cert)
        % set_certainty(cert) - Sets the certainty of computations involving
        % parameters x, y, and lambda for which we have std deviations.
        % assuming all errors are normally distributed we use the std_dev
        % to adjust all bounds to be accurated within probaiblity=certainty
        % (in) cert - 0<float<1: The certainty with wich to compute all 
        %                        error depended calculations. (default=0.95)
        assert(0<cert && cert<1);
        obj.certainty=cert;
        obj.normal_deviation=norminv(cert);
    end

    function lengths=trackLengths(obj, idxs)
        % lengths=obj.trackLengths() - Return the lengths (number of localizations) of each
        % track.
        % (in) idxs (optional) - An index array or logical array of track indexes to consider.  
        %      Default: consider all tracks.
        % (out) lengths - A nTx1 vector of track lengths
        if nargin==1
            Ts=obj.T;
        elseif nargin==2
            Ts=obj.T(idxs);
        end
        lengths=cellfun(@(t) size(t,1), Ts);
    end
    
    function times=trackTimeSpans(obj, idxs)
        % times=obj.trackTimeSpans() - Return the start and end time in (s) for each track.
        % (in) idxs (optional) - An index array or logical array of track indexes to consider.  
        %      Default: consider all tracks.
        % (out) times - A nTx2 matrix where col 1 is track start time and
        % col 2 is track end time.
        if nargin==1
            Ts=obj.T;
        elseif nargin==2
            Ts=obj.T(idxs);
        end
        times=cellmatfun(@(t) [t(1,1) t(end,1)], Ts);
    end

    function durations=trackDurations(obj, idxs)
        % durations=obj.trackDurations() - Returns the duration (s) of each track.
        % (in) idxs (optional) - An index array or logical array of track indexes to consider.  
        %      Default: consider all tracks.
        % (out) durations - A nTx1 vector of durations
        if nargin==1
            ts=obj.trackTimeSpans();
        elseif nargin==2
            ts=obj.trackTimeSpans(idxs);
        end
        durations=ts(:,2)-ts(:,1);
    end

    function bounds=trackSpectraBounds(obj, idxs)
        % bounds=obj.trackSpectraBounds() - Return the maximum range of
        % the track wavelength assuming the errors are normally distributed and the
        % std_lambda values are given.
        % (in) idxs (optional) - An index array or logical array of track indexes to consider.  
        %      Default: consider all tracks.
        % (out) bounds - a nTx2 matrix where col 1 is the lower-bound and col 2 the upper-bound
        if nargin==1
            Ts=obj.T;
        elseif nargin==2
            Ts=obj.T(idxs);
        end
        bounds=cellmatfun(@obj.computeSpectraBound,Ts);
    end

    function bboxes=trackBBoxes(obj, idxs)
        % bbox=obj.trackBBoxes() - Return the maximum x and y ranges of the
        % track localization positions assuming the errors are normally distributed
        % and the std_x and std_y are given.
        % (in) idxs (optional) - An index array or logical array of track indexes to consider.  
        %      Default: consider all tracks.
        % (out) bboxs - a nTx4 matrix where: col1:xmin col2:xmax col3:ymin col4:ymax
        if nargin==1
            Ts=obj.T;
        elseif nargin==2
            Ts=obj.T(idxs);
        end
        bboxes=cellmatfun(@obj.computeBBox,Ts);
    end

    %%% Group Statistics Methods %%%

    function time=groupTimeSpan(obj)
        % time=obj.groupTimeSpan() - Return the earliest start and latest end time in (s) over all tracks
        % in this group.
        % (out) time - A 1x2 matrix where col 1 is earliest start time and
        % col 2 is latest end time.
        L=obj.Localizations;
        time=[min(L(:,1)) max(L(:,1))];
    end
    
    function duration=groupDuration(obj)
        % duration=obj.groupDuration() - Returns the overall duration (s) of all tracks in the group.
        % This is the timespan from the start of the earliest track to the end of the latest track
        % (out) duration - A scalar giving the duration
        tm=obj.groupTimeSpan();
        duration=tm(2)-tm(1);
    end

    function bbox=groupBBox(obj)
        % bbox=obj.groupBBox() - Return the maximum x and y ranges of all tracks
        % localization positions from this group assuming the errors are normally distributed
        % and the std_x and std_y are given.
        % (out) bbox - a 1x4 matrix where: col1:xmin col2:ymin col3:xmax col4:ymax
        bbox=obj.computeBBox(obj.Localizations);
    end
    
    function bound=groupSpectraBound(obj)
        % bound=obj.groupSpectraBound() - Return the maximum range of
        % the wavelength range over all tracks assuming the errors are normally distributed and the
        % std_lambda values are given.
        % (out) bounds - a 1x2 vector where col 1 is the lower-bound of lambda over all tracks
        % and col 2 the upper-bound of lambda over all tracks
        bound=obj.computeSpectraBound(obj.Localizations);
    end

    %%% Track Filtering Methods %%%

    function nrejected=filterLength(obj, min_length)
        % nrejected=obj.filterLength(min_length) - Reject tracks with fewer than min_length Localizations.
        % All tracks not meeting restriction are removed along with their localizations.
        % (in) min_length - int>0: minimum length of a track (in localization count)    
        % (out) nrejected - int>=0: number of rejected tracks
        length=obj.trackLengths();
        nrejected=obj.filterSelection( length<min_length );
    end

    function nrejected=filterDuration(obj, min_duration)
        % nrejected=obj.filterDuration(min_duration) - Reject tracks with duration less than min_duration
        % All tracks not meeting restriction are removed along with their localizations.
        % (in) min_duration - float>0: minimum duration (s) of a track.
        % (out) nrejected - int>=0: number of rejected tracks
        duration=obj.trackDurations();
        nrejected=obj.filterSelection( duration<min_duration );
    end

    function nrejected=filterSpectralDeviation(obj, max_spectral_deviation)
        % nrejected=obj.filterSpectralDeviation(max_spectral_deviation) - Reject tracks with
        % spectral deviation greater than max_spectral_deviation with given certainty, assuming errors are
        % normally distributred with and track std_lambda is accurate.
        % All tracks not meeting restriction are removed along with their localizations.
        % (in) max_spectral_deviation - float>0: maximum spectral width of track in (micron).
        % (out) nrejected - int>=0: number of rejected tracks
        bounds=obj.trackSpectraBounds();
        deviation=bounds(:,2)-bounds(:,1);
        nrejected=obj.filterSelection( deviation>max_spectral_deviation );
    end

    function nrejected=filterRandomSubsample(obj, keep_fraction)
        [~,todrop]=sort(rand(obj.nT, 1));
        n=round(keep_fraction*obj.nT)+1;
        nrejected=obj.filterSelection(todrop(n:end));
    end
   
    function nrejected=filterSelection(obj, idxs)
        % nrejected=obj.filterSelection(idxs) - Reject tracks with given indexs.
        % All tracks in idxs are removed along with their localizations.  This provides a way to remove
        % custom subsets of tracks based on custom criteria.
        % (in) idxs - index vector or logical vector of tracks to remove.  All indexs must be valid and no repeats!
        % (out) nrejected - int>=0: number of rejected tracks
        nrejected=obj.nT;
        obj.T(idxs)=[];
        nrejected=nrejected-obj.nT;            
    end
    


    function track_table=trackTable(obj, idxs)
        if nargin==1
            Ts=obj.T;
        elseif nargin==2
            Ts=obj.T(idxs);
        end
        track_table=cellmap(@obj.makeTrackTable,Ts);
    end

    
    function f=plot_localizations(obj,idxs)
        if nargin<2
            idxs=1:obj.nLocalizations;
        end
        L=obj.Localizations(idxs,:);
        f=figure('Renderer','OpenGL');
        rgb=wavelengthRGB(L(:,4), obj.groupSpectraBound());
        ms=10;
        scatter3(L(:,2), L(:,3), L(:,1),ms,rgb,'fill');
        title('Localizations Colored by Lambda');
        zlabel('Time (frames)');
        xlabel('X (pixels)');
        ylabel('Y (pixels)');
    end

    function f=plot_tracks(obj,idxs)
        if nargin<2
            idxs=1:obj.nT;
        end
        f=figure;
        ms=10;
        hold;
        for i=idxs
            t=obj.T{i};
            rgb=wavelengthRGB(t(:,4),obj.groupSpectraBound());
            scatter3(t(:,2),t(:,3),t(:,1),ms,rgb,'fill');
            plot3(t(:,2), t(:,3), t(:,1),'k-');
        end
        hold;
        title('Tracks')
        zlabel('Time (frames)');
        xlabel('X (pixels)');
        ylabel('Y (pixels)');
    end
end %public methods

properties (Access=private)
    %This is calculated based on the obj.certainty property.
    % it is norminv(certainty)
    normal_deviation;
end %private properties

methods (Access=private)
    function bbox=computeBBox(obj,track)
        % bbox=obj.computeBBox(track) - compute a bbox that contains all x,y positions of
        % the particle at all time with high probability.  The bbox is always aligned with
        % the xy axes, and spans all localizations for the particle.
        % (in) track - A track from obj.T
        % (out) bbox - 1x4 double: [minx, miny, maxx, maxy] /w obj.certainty
        minx=min(track(:,2)-track(:,6)*obj.normal_deviation);
        miny=min(track(:,3)-track(:,7)*obj.normal_deviation);
        maxx=max(track(:,2)+track(:,6)*obj.normal_deviation);
        maxy=max(track(:,3)+track(:,7)*obj.normal_deviation);
        bbox=[minx miny maxx maxy];
    end
    
    function bound=computeSpectraBound(obj,track)
        % bbox=obj.computeBBox(track) - compute a range that contains all lambda mean values
        % with with high probability over all localizations for the particle.
        % (in) track - A track from obj.T
        % (out) bound - 1x2 double: [min_lambda, max_lambda] /w obj.certainty
        minval=min(track(:,4)-track(:,8)*obj.normal_deviation);
        maxval=max(track(:,4)+track(:,8)*obj.normal_deviation);
        bound=[minval maxval];
    end

    function track_table=makeTrackTable(obj,track)
        % track_table=makeTrackTable(track) - make a table from a track matrix.  This uses
        % the constant member variables TColumns, TUnits, etc. to make a fancy (but slow)
        % table object.
        % (in) track - A track from obj.T
        % (out) - A fancy table!
        track_table=array2table(track);
        track_table.Properties.VariableNames=obj.TColumns;
        track_table.Properties.VariableUnits=obj.TUnit;
        track_table.Properties.VariableDescriptions=obj.TDescription;
        track_table.Properties.Description=obj.TTitle;
        track_table.Properties.DimensionNames=obj.TDimensionNames;
    end

    function initialize_SPTHSI(obj, spthsi)
        % obj.initialize_SPTHSI(spthsi) - Given an spthsi object, intializes the track and localization data structures
        %
        % Tracks are ordered by length. And all units are converted to nm or s.
        obj.filename=spthsi.DataFile;
        fparts=strsplit(obj.filename,{'\\','/'});
        [~,obj.name,~]=fileparts(fparts{end});
        nT=numel(spthsi.Tracks);
        obj.T=cell(nT,1);
        lengths=arrayfun(@(t) numel(t.X), spthsi.Tracks); %Get length of each track
        fprintf('Total Tracks: %i\n',nT);
        fprintf('Total Localizations: %i\n',sum(lengths));

        %extract pixel and frame units from spthsi
        psize=spthsi.Stats.Data.acqParams.pixelSize*1e-3; %This should be pixel size in nm (convert to microns)
        ftime=spthsi.Stats.Data.acqParams.t; %This should be frame aquisition time in s.
        obj.frameT=ftime;

        %copy tracks to obj.T
        [~,sidx]=sort(lengths,'descend'); %get indexs of sorted tracks
        for i = 1:nT
            ts=spthsi.Tracks(sidx(i));
            obj.T{i}=[ftime*(ts.Frame'-1), psize*ts.X', psize*ts.Y', ts.Wv', ts.Photons',...
                      psize*psize*(ts.std_x.*ts.std_x)', psize*psize*(ts.std_y.*ts.std_y)', (ts.std_w.*ts.std_wv)',...
                      (ts.std_Photons.*ts.std_Photons)' ts.Frame'];
        end
    end

end %private methods
methods (Static=true)
   function [ Xp ] = AffineTform(X,M,Xref,NN)
        %AffineTform Affine Transform of a set of points
        %   AffineTform(X,M) performs the Affine transform of X to Xp
        %   by M*X=Xp where the M matrix has the structure:
        %
        %   [a b c][x] [xp]
        %   [d e f]|y]=[yp]
        %   [0 0 1][1] [1 ]
        %
        %   X must be Nx2
        %
        %   AffineTform(X,M,Xref) performs a local transform using a 1/distance
        %   weighted average of the NN closest control points.
        %   M must be 3x3xN and Xref must be Nx2.

        if nargin<3
            [Xp]=AFT(X,M);
            return
        end

        Npoints=size(X,1);
        Xp=zeros(size(X));

        for ii=1:Npoints
            r2=(Xref(:,1)-X(ii,1)).^2+(Xref(:,2)-X(ii,2)).^2;
            r=sqrt(r2);
            [val idx]=sort(r2);
            x=0;
            if r(idx(1))==0
                [x]=AFT(X(ii,:),M(:,:,idx(1))); 
            else
                r0=sum(1./r(idx(1:NN)));
                for nn=1:NN
                    id=idx(nn);
                    [xtmp]=AFT(X(ii,:),M(:,:,id));
                    w=1/r(id)/r0;
                    x=x+xtmp*w;
                end

            end

            Xp(ii,:)=x;
        end
        function [Xp]=AFT(X,M)
            xtmp=cat(2,X,ones(size(X,1),1));
            xfix =  xtmp*M(1,:)';
            yfix =  xtmp*M(2,:)';
            Xp=[xfix yfix];
        end
   end
   function [M Minv] = FindAffineTform(X,Xp,NN)
        %FindAffineTform Find Affine Transformation Matrix
        %   FindAffineTform(X,Xp) finds the Affine transform of X to Xp
        %   by M*X=Xp where the M matrix has the structure:
        %
        %   [a b c][x] [xp]
        %   [d e f]|y]=[yp]
        %   [0 0 1][1] [1 ]
        %
        %   X and Xp must be Nx2
        %
        %   FindAffineTform(X,Xp,NN) uses the closest NN points to each elelement
        %   in Xp to perform a local transformation, giving 3x3xN output matrix

        if nargin<3
            [M Minv]=FATM(X,Xp);
            return
        end

        Npoints=size(Xp,2);
        ATM=zeros(3,3,Npoints);
        ATMinv=zeros(3,3,Npoints);

        for ii=1:size(Xp,1)
            r2=(Xp(:,1)-Xp(ii,1)).^2+(Xp(:,2)-Xp(ii,2)).^2;
            [val idx]=sort(r2);
            XpNN=Xp(idx(1:NN),:);
            XNN=X(idx(1:NN),:);
            [Mtmp Minvtmp] = FindAffineTform(XNN,XpNN);
            M(:,:,ii)=Mtmp;
            Minv(:,:,ii)=Minvtmp;
        end



        function [M Minv]=FATM(X,Xp)

            A=cat(2,X,ones(size(X,1),1));
            invA=pinv(A);
            ATmatrixRow1=invA*Xp(:,1);
            ATmatrixRow2=invA*Xp(:,2);

            M=cat(1,ATmatrixRow1',ATmatrixRow2',[0 0 1]);
            Minv=pinv(M);

        end
    end

end %Static Methods


end %class

