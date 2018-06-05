classdef PairAnalysis < Pickle & GUIBuilder

    % PairAnalysis - Represents pairs of tracks, which come from a track group.  A
    % PairAnalysis object can be formed from multiple TrackGroups which each represent
    % seperate (or subsequent) experiments.  A individual pair my only come from a single track
    % group.
    %
    % Design decisions:
    % * As with TrackGroup, I wish we could make P a cell array of tables, but it is at present over 1000
    %   times slower than a cell array of matricies.  If parameters are added to the pair
    %   matrix, it will not be pleasent to updates everything cleanly.  And for 3D data or
    %   other non-standard tracking, a new class or subclass will be needed.
    %
    %
    % Author: Mark J. Olah (mjo@cs.unm.edu)
    % Date: Feb-2014
    %

    properties (Constant=true, Hidden=true)



        % The primary purpose of the PairAnalysis base class is to provide a common format for paired-track data.
        % The data for each pair can be retrieved in multiple formats simmilar to the useage of the RPT.
        % These constant properties describe the data format, where each pair 

        NPairColumns = 18;  %Number of columns in a pair matrix

        % Column names of a pair matrix (or table)
        %       (L=lambda) is the wavelength in nm.
        PairColumnNames = {'t',... % 1
                           'x1', 'y1', 'L1', 'I1',... % 2 3 4 5 
                           'x2', 'y2', 'L2', 'I2',... % 6 7 8 9
                           'SE_x1', 'SE_y1', 'SE_L1', 'SE_I1',... % 10 11 12 13
                           'SE_x2', 'SE_y2', 'SE_L2', 'SE_I2',... % 14 15 16 17
                           'frameIdx'};% 18

        % Units of pair table columns
        PairColumnUnits = {'s',... % 1
                           'um', 'um', 'nm', '#',... % 2 3 4 5 
                           'um', 'um', 'nm', '#',... % 6 7 8 9
                           'um', 'um', 'nm', '#',... % 10 11 12 13
                           'um', 'um', 'nm', '#',... % 14 15 16 17
                           '#'}; %18

        % Description of pair table columns
        PairColumnDescriptions = {'time',... % 1
                                  'p1 x-position', 'p1 y-position', 'p1 mean-lambda', 'p1 intensity',... % 2 3 4 5 
                                  'p2 x-position', 'p2 y-position', 'p2 mean-lambda', 'p2 intensity',... % 6 7 8 9
                                  'p1 x standard deviation', 'p1 y standard deviation', 'p1 lambda standard deviation', 'p1 intensity standard deviation',... % 10 11 12 13
                                  'p2 x standard deviation', 'p2 y standard deviation', 'p2 lambda standard deviation', 'p2 intensity standard deviation',... % 14 15 16 17
                                  'frame index(1-based)'};% 18

        PairTableTitle = 'PairLocalizations';  %Title of a pair table   
        PairTableDimensionNames = {'Localization', 'Coordinates'}; %DimensionNames of a pair table 

        %The distance matrix/table stores computed distances (mean distance min distance and max distance) for each
        %point in each trajectory-pair.  This allows rapid plotting and retrieval of the distance information.
        NDistColumns = 13; % Number of columns in a distance matrix
        DistColumnNames = {'t',... % 1 
                           'd', 'd_min', 'd_max',... % 2 3 4
                           'dx', 'dx_min', 'dx_max',... % 5 6 7
                           'dy', 'dy_min', 'dy_max',... % 8 9 10
                           'dL', 'dL_min','dL_max'} % 11 12 13

        DistColumnUnits = {'s',... % 1
                           'um', 'um', 'um',... % 2 3 4
                           'um', 'um', 'um',... % 5 6 7
                           'um', 'um', 'um',... % 8 9 10
                           'nm', 'nm', 'nm'} % 11 12 13

        DistColumnDescriptions = {'time',... % 1
                                  '2D dist', '2D dist lower bound', '2D dist upper bound',... % 2 3 4
                                  'x-dist', 'x-dist lower bound', 'x-dist upper bound',... % 5 6 7
                                  'y-dist', 'y-dist lower bound', 'y-dist upper bound',... % 8 9 10
                                  'lambda-dist', 'lambda dist lower bound', 'lambda dist upper bound'} % 11 12 13
        DistTableTitle = 'PairDistances';
        DistTableDimensionNames = {'Localization', 'Distances'};
    end % public constant properties

    properties
        
        physicalBounds; % [xmin, xmax, ymin, ymax] in micron.  The smallest axis-aligned rectangular bounding box than encompasses bothe of the ChannelPhysicalBBox's
        timeBounds; % [tmin, tmax] in seconds.
        frameBounds; % [frame min frame max]
        frameT; % Frame time, retrieved from obj.data
        

        tracks; % A cell-array of tracks.  The format is determined by the subclass as HS and 2-color data will be stored differently

        %pairsIdentified: Boolean true if potential pairs have been identified already.  
        % This is set by calling identifyPairs after tracks have been filtered appropriately.  
        % Pairs can then be further filtered using the filterPair...() methods
        pairsIdentified = false;
        
        % pairs: A cell-array of pair matricies.  Each matrix has a row for each common localization
        pairs;
        
        % pairIds: nPairsx2 matrix of indexes.  Each row represents a pair giving the indexes into the tracks
        % cell-array for each of the two tracks in the pair.
        pairIds;
        
        % pairStats: struct array of computed statistics for each pair in the obj.pairs datastructure
        pairStats;
        
        % pairDists: cell-array size:[nPairs,1]. Distance statisitcs for each pair
        pairDists;
    end %public properties

    properties (Abstract = true)
        minPairLength; %Minimum lenght for a pair to be considered valid
        minPairDistance; %Minimum approach distance necessary for a pair to be considered valid
    end

    properties (SetAccess = protected, GetAccess = public)
        pair_min_length;   % The minimum pair length (# of simultaneous localizations) uses in pair idetification
        pair_min_distance; % The minumum approach distance used as a cutoff in the pair identification
        certainty=0.95;     % The probability threshold to which comparisons and calculations involving uncertainty are performed.
    end

    properties (Transient=true, Hidden=true)
        normal_deviation;
    end

    properties (Transient=true)
        data;       %Link to the SPData or HSData corresponding to these tracks
    end   

    properties (Dependent=true)      
        nPairs;  % number of potential pairs identified
        nTracks; % number of tracks
    end 

    properties (Hidden=true)
        version=1; %For future file format version changes
    end

methods
    function obj=PairAnalysis()
        % PairAnalysis() - Make an empty pair analysis object
        % PairAnalysis(name) -
        %    (in) name - a name for this PairAnalysis instance
        if nargin==1
            obj.name=name;
        end
    end
    
    function clearPairs(obj)
        % Clear all pairs and related pair data structures.  Sets obj.pairsIdentified=false
        obj.pairs={};
        obj.pairIds={};
        obj.pairStats={};
        obj.pairDists={};
        obj.pairsIdentified=false;
    end

    function pBounds = pairPhysicalBounds(obj, Ps)
        % [in]
        %  Ps - A cell array of pari matricies in PairAnalysis format.
        % [out]
        %  pBounds = [minx, maxy, miny, maxy] (microns)
        Ps = makecell(Ps);
        D = obj.normal_deviation;
        minx = min(cellfun(@(P) min([P(:,2)-P(:,10)*D; P(:,6)-P(:,14)*D]), Ps));
        miny = min(cellfun(@(P) min([P(:,3)-P(:,11)*D; P(:,7)-P(:,15)*D]), Ps));
        maxx = max(cellfun(@(P) max([P(:,2)+P(:,10)*D; P(:,6)+P(:,14)*D]), Ps));
        maxy = max(cellfun(@(P) max([P(:,3)+P(:,11)*D; P(:,7)+P(:,15)*D]), Ps));
        pBounds = [minx maxx miny maxy];
    end  


end %public methods

methods % Dependent property accessor methods
    function nP = get.nPairs(obj)
        nP = numel(obj.pairs);
    end
    
    function nT = get.nTracks(obj)
            nT = numel(obj.tracks);
    end

    function dev = get.normal_deviation(obj)
        % The number of sigmas out we should take for error bars to give obj.certainty chance of
        % containing true value.  This assumes a guassian error distribution.
        dev = norminv(obj.certainty);
    end
end %dependent properties

methods (Access=protected, Abstract=true)
    Cs = setTrackColors(obj, Tidxs, method);
end

methods (Access=protected)
    function Tidxs = getPairTrackIdx(obj, Pidxs)
        Tidxs = union(obj.pairIds(Pidxs,1), obj.pairIds(Pidxs,2));
    end

    function windows = computeProximityWindows(obj, pairDist, dist_threshold)
        dist_ok_delta = diff([false; pairDist(:,4)<=dist_threshold; false]);
        start_times = dist_ok_delta==1;
        end_times = dist_ok_delta==-1;
        times = [pairDist(:,1); pairDist(end,1)+obj.frameT];
        windows = [times(start_times) times(end_times)];
    end

    function ptable = makePTable(obj, pair)
        ptable = array2table(pair);
        ptable.Properties.VariableNames = obj.PColumn;
        ptable.Properties.VariableUnits = obj.PUnit;
        ptable.Properties.VariableDescriptions = obj.PDescription;
        ptable.Properties.Description = obj.PTitle;
        ptable.Properties.DimensionNames = obj.PDimensionNames;
    end

    function dtable = makeDTable(obj, dist)
        dtable = array2table(dist);
        dtable.Properties.VariableNames = obj.DColumn;
        dtable.Properties.VariableUnits = obj.DUnit;
        dtable.Properties.VariableDescriptions = obj.DDescription;
        dtable.Properties.Description = obj.DTitle;
        dtable.Properties.DimensionNames = obj.DDimensionNames;
    end
end % private methods

methods (Static=true)
    function checkPairFormat(pair)
        % Check that the pair matrix is in valid format
        nObs = size(pair,1);
        if nObs<1
            error('PairAnalysis:checkPairFormat','Pair has no observations');
        end
        if size(pair,2) ~= PairAnalysis.NPairColumns 
            error('PairAnalysis:checkPairFormat','Pair has %i columns.  Expected %i.',size(pair,2),PairAnalysis.NPairColumns);
        end
        times = pair(:,1);
        if any(times<0 | ~isfinite(times))
             error('PairAnalysis:checkPairFormat','Times must be real values >0');
        end
        if any(sort(times) ~= times)
             error('PairAnalysis:checkPairFormat','Times must be increasing');
        end
        if numel(unique(times)) ~= numel(times)
             error('PairAnalysis:checkPairFormat','Times must be unique');
        end
    end

    function stepcor = computeStepAngularCorrelation(pair)
        s1 = pair(2:end,[2,3])-pair(1:end-1,[2,3]);
        s2 = pair(2:end,[6,7])-pair(1:end-1,[6,7]);
        stepcor = dot(s1,s2,2)./sqrt(dot(s1,s1,2).*dot(s2,s2,2));
    end

    function dists=makePairDists(pair, certainty)
        if nargin==1
            certainty=0.95;
        end
        normal_deviation=norminv(certainty);
        dists=zeros(size(pair,1), PairAnalysis.NDistColumns);
        
        %dx
        dist_x=abs(pair(:,2)-pair(:,6));
        delta_x=(pair(:,10)+pair(:,14)).*normal_deviation;

        %dy
        dist_y=abs(pair(:,3)-pair(:,7));
        delta_y=(pair(:,11)+pair(:,15)).*normal_deviation;

        %d lambda
        dist_lambda=abs(pair(:,4)-pair(:,8));
        delta_lambda=(pair(:,12)+pair(:,16)).*normal_deviation;
        
        %d
        dist_2D=sqrt(dist_x.^2+dist_y.^2);
        delta_2D= (sqrt( dist_x.^2 .* (pair(:,10).^2+pair(:,14).^2) + ...
                         dist_y.^2 .* (pair(:,11).^2+pair(:,15).^2))./dist_2D).*normal_deviation;

        %fill in dists array
        dists(:,1)=pair(:,1);
        dists(:,2)=dist_2D;
        dists(:,3)=max(0, dist_2D - delta_2D);
        dists(:,4)=dist_2D + delta_2D;
        dists(:,5)=dist_x;
        dists(:,6)=max(0, dist_x - delta_x);
        dists(:,7)=dist_x + delta_x;
        dists(:,8)=dist_y;
        dists(:,9)=max(0, dist_y - delta_y);
        dists(:,10)=dist_y + delta_y;        
        dists(:,11)=dist_lambda;
        dists(:,12)=max(0, dist_lambda - delta_lambda);
        dists(:,13)=dist_lambda + delta_lambda;
    end 
end


end %classdef

