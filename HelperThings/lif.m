classdef lif < handle
    %  LIF
    %  class for loading data from lif (leica image format) files
    %  not sure whether it works for z-stacks, didn't have test data
    %  
    %  SYNTAX:
    %  lifObj = lif(fileDir,fileName)
    %
    %  NEEDED:
    %  Bio-Formats toolbox for Matlab, download from:
    %  http://downloads.openmicroscopy.org/bio-formats/5.1.3/
    %
    %  Marjolein Meddens 2015 UNM
    
    properties
        FileDir % directory of lif file
        FileName % name of lif file
    end
    
    properties (SetAccess=protected)
        NumSeries % number of series in the lif file
        SeriesName % cell array of series names
        SeriesSize % [x, y, z, t, c]
    end
    
    properties (SetAccess=protected, Hidden=true)
        Reader
    end
    
    methods
        
        function obj = lif(FileDir,FileName)
            % constructor
            if isempty(strfind(FileName,'.lif'))
                FileName = [FileName '.lif'];
            end
            if ~exist(fullfile(FileDir,FileName),'file')
                error('lif:FileNotFound','lif file not found');
            else
                obj.FileDir = FileDir;
                obj.FileName = FileName;
                obj.Reader = bfGetReader(fullfile(FileDir,FileName));
            end
            obj.getSeriesInfo;
        end
        
        function getSeriesInfo(obj)
            % LIF.GETSERIESINFO
            % Method to get name and file size of all series in lif file
            % Will populate NumSeries, SeriesName and SeriesSize properties
            
            obj.NumSeries = obj.Reader.getSeriesCount();
            omeMeta = obj.Reader.getMetadataStore(); % contains ome metadata for all series
            % preallocate variables
            Series = cell(obj.NumSeries,1);
            Sizes = zeros(obj.NumSeries,5);
            % loop through series
            for ii = 1 : obj.NumSeries
                obj.Reader.setSeries(ii-1); % set series to be able to get lif metadata for current series
                lifMetaData = obj.Reader.getSeriesMetadata();
                Series{ii} = lifMetaData.get('Image name');
                szX = omeMeta.getPixelsSizeX(ii-1).getValue(); % image width, pixels
                szY = omeMeta.getPixelsSizeY(ii-1).getValue(); % image height, pixels
                szZ = omeMeta.getPixelsSizeZ(ii-1).getValue(); % number of Z slices
                szT = omeMeta.getPixelsSizeT(ii-1).getValue(); % number of Time points
                szC = omeMeta.getPixelsSizeC(ii-1).getValue(); % number of Channels
                Sizes(ii,:) = [szX, szY, szZ, szT, szC];
            end
            obj.SeriesName = Series;
            obj.SeriesSize = Sizes;
        end
        
        function [sequence, metadata] = loadSeries(obj,seriesNum)
            % LIF.LOADSERIES(seriesNum)
            % Method to load a specific series from lif file
            %
            % INPUT
            % seriesIdx:    index of series to be loaded
            %               corresponds to order in lifObj.SeriesName
            % 
            % OUTPUT
            % sequence:     5D array containing image data [x y z t c]
            %               same data type as original data
            % metadata:     structure containing series metadata
            %               some specific info has been subtracted from lif
            %               all metadata is found as string in omeMetaData
            %               and 
            %
            % NOTE:
            % let me know if you want other specific metadata to be
            % retrieved, it's easy to get more
            % Marjolein
            
            % get metadata from the reader
            obj.Reader.setSeries(seriesNum-1); % set series
            % ome metadata
            omeMeta = obj.Reader.getMetadataStore(); % contains ome metadata for all series
            pxSzX = omeMeta.getPixelsPhysicalSizeX(seriesNum-1).value(ome.units.UNITS.MICROM);
            metadata.PxSzX = pxSzX.doubleValue(); % um
            pxSzY = omeMeta.getPixelsPhysicalSizeY(seriesNum-1).value(ome.units.UNITS.MICROM);
            metadata.PxSzY = pxSzY.doubleValue(); % um
            if obj.SeriesSize(seriesNum,3)>1
                pxSzZ = omeMeta.getPixelsPhysicalSizeZ(seriesNum-1).value(ome.units.UNITS.MICROM);
                metadata.PxSzZ = pxSzZ.doubleValue(); % um
            end
            pixType = omeMeta.getPixelsType(seriesNum-1);
            metadata.dataType = char(pixType.getValue);
            % leica metadata
            leicaMetaData = obj.Reader.getSeriesMetadata();
            metadata.frameTime = str2double(leicaMetaData.get('ATLConfocalSettingDefinition|FrameTime')); % sec
            metadata.objectiveName = leicaMetaData.get('ATLConfocalSettingDefinition|ObjectiveName');
            metadata.pixelDwellTime = str2double(leicaMetaData.get('ATLConfocalSettingDefinition|PixelDwellTime')); % sec
            % dump all metadata as strings
            metadata.omeMetaData = char(omeMeta.dumpXML()); % string of all ome metadata
            metadata.leicaMetaData = char(leicaMetaData); % string of all leica metadata
            % load the sequence
            % preallocate sequence
            sequence = zeros(obj.SeriesSize(seriesNum,:),metadata.dataType);
            % retrieve all planes
            for cc = 1 : obj.SeriesSize(seriesNum,5)
                for tt = 1 : obj.SeriesSize(seriesNum,4)
                    for zz = 1 : obj.SeriesSize(seriesNum,3)
                        iPlane = obj.Reader.getIndex(zz-1, cc-1, tt-1) + 1;
                        sequence(:,:,zz,tt,cc) = bfGetPlane(obj.Reader, iPlane);
                    end
                end
            end
        end
        
    end
        
end