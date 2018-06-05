% File: Tracker.m
%
% Mark J. Olah (mjo@cs.unm.edu)
% 04/2015

classdef Tracker < IfaceMixin
    properties (Constant=true)
        Trackers = {'LAPTrack2D'};
    end

    properties (SetAccess=protected)
        trackerType;
    end

    properties (Access = protected, Transient=true)
        initialized=false; %True once the object is correctly initialized by the initializeTrack method
    end
    
    methods
        function obj=Tracker(trackerType,params)
            %
            iface = str2func(sprintf('%s_Iface',trackerType));
            obj=obj@IfaceMixin(iface);
            obj.trackerType = trackerType;
            obj.initialized=obj.openIface(params);
        end
        
        function stats = getStats(obj)
            stats=obj.call('getStats');
        end

        function initializeTracks(obj, frameIdx, position, positionSE, feature, featureSE)
            frameIdx = int32(frameIdx(:));
            if nargin==4
                obj.call('initializeTracks',frameIdx, position, positionSE);
            elseif nargin==6
                obj.call('initializeTracks',frameIdx, position, positionSE, feature, featureSE);
            end
        end

        function [curIdx, nextIdx, costMat, connections, conn_costs] = debugF2F(obj, frameIdx)
            [curIdx, nextIdx, costMat, connections, conn_costs] = obj.call('debugF2F',int32(frameIdx));
        end

        function [costMat, connections, conn_costs] = debugCloseGaps(obj)
            [costMat, connections, conn_costs]= obj.call('debugCloseGaps');
        end
        
        function nTracks = linkF2F(obj)
            nTracks = obj.call('linkF2F');
        end
        
        function nTracks = closeGaps(obj)
            nTracks = obj.call('closeGaps');
        end

        function tracks = generateTracks(obj)
            tracks = obj.call('generateTracks');
        end
        
        function tracks = getTracks(obj)
            tracks = obj.call('getTracks');
        end
    end %public methods
    methods (Access=protected)
        function checkPoints(obj,points)
            if size(points,2) ~=5
                error('SRRender2D:checkPoints','Points Incorrect number of columns');
            end
            if any( points(:,1)<=0 )
                error('SRRender2D:checkPoints','Bad intensity values');
            end
            if any( points(:,2)<0 | points(:,2)>obj.size(1) )
                error('SRRender2D:checkPoints','Bad X values');
            end
            if any( points(:,3)<0 | points(:,3)>obj.size(2) )
                error('SRRender2D:checkPoints','Bad Y values');
            end
        end
    end %protected methods
end %classdef
