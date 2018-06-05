% File: KDTree.m
%
% Mark J. Olah (mjo@cs.unm.edu)
% 09/2014

classdef KDTree < IfaceMixin
    properties (SetAccess=protected)
        Ndim; %Number of dimensions
        N;    %Number of points
        datatype;
        datacaster; %Handle to method to convert data to this type
        bbox; %The bounding box for all data.  size:=[Ndim 2]
              %  col 1 minimum corner for each dimension 
              %  col 2 maximum corner for each dimension 
    end
    
    properties (Access = protected, Transient=true)
        initialized=false; %True once the object is correctly initialized by the initializeTrack method
    end
    
    methods
        function obj=KDTree(points)
            % Construct a new KDTree object or a blank object.
            %
            % Usage:
            %   obj=KDTree(points);
            %
            % Inputs:
            %   points=[Ndim X N] matrix of Ndim dimensional points.  With N total points.
            %
            if isa(points,'double')
                iface=@KDTreeDouble_Iface;
                datatype='double';
                datacaster=@double;
            elseif isa(points,'single')
                iface=@KDTreeSingle_Iface;
                datatype='single';
                datacaster=@single;
            elseif isa(points,'int32')
                iface=@KDTreeInt_Iface;
                datacaster=@int32;
                datatype='int32';
            else
                error('KDTree:constructor','Unsupported type: %s', class(a));
            end
            obj=obj@IfaceMixin(iface);
            obj.datatype=datatype;
            obj.datacaster=datacaster;
            obj.Ndim=size(points,1);
            obj.N=size(points,2);
            assert(obj.N>0 && obj.Ndim>0);
            obj.bbox=[min(points,[],2), max(points,[],2)];
            obj.initialized=obj.openIface(points);
        end 
        
        function initialize(obj,points)
            % Inputs:
            %   points=[Ndim X N] matrix of Ndim dimensional points.  With N total points.  
            %                    Must be same type as obj.datatype
            %
            if obj.initialized
                obj.closeIface();
            end
            obj.initialized=false;
            if ~strcmp(class(points),obj.datatype)
                error('KDTree:initialize','Cannot change type of points in reinitialization.');
            end
            obj.Ndim=size(points,1);
            obj.N=size(points,2);
            assert(obj.N>0 && obj.Ndim>0);
            obj.bbox=[min(points,[],2), max(points,[],2)];
            obj.initialized=obj.openIface(points);
        end
        
        function points=query(obj, min_corner, max_corner)
            %
            %
            % Useage:
            %    points=obj.query(min_corner, max_corner);
            %
            % Inputs:
            %  min_corner - size:[Ndim] Location of query range corner with minimum coordinates
            %  max_corner - size:[Ndim] Location of query range corner with maximum coordinates
            %
            % Ouptuts:
            %  points - [Ndim, M] Where M is the number of points in the query range
            %
            if ~isvector(min_corner) || length(min_corner)~=obj.Ndim
                error('KDTree:query', 'min_corner should be vector of length Ndim=%i',obj.Ndim);
            end
            if ~isvector(max_corner) || length(max_corner)~=obj.Ndim
                error('KDTree:query', 'max_corner should be vector of length Ndim=%i',obj.Ndim);
            end
            points=obj.call('query', obj.datacaster(min_corner), obj.datacaster(max_corner));
        end
        
        function points=queryCircle(obj, center, radius)
            % Find all points in given circle
            %
            % Useage:
            %    points=obj.queryCircle(center, radius);
            %
            % Inputs:
            %  center - size:[Ndim] Center of circle 
            %  radius - scalar>0 Radius of circle
            %
            % Ouptuts:
            %  points - [Ndim, M] Where M is the number of points in the query circle
            %
            if ~isvector(center) || length(center)~=obj.Ndim
                error('KDTree:queryCircle', 'Center is must be a vector of length Ndim=%i',obj.Ndim);
            end
            if radius<=0
                error('KDTree:queryCircle', 'Radius must be positive');
            end
            points=obj.call('query', obj.datacaster(center-radius),  obj.datacaster(center+radius));
            points=points(:, sum((points-repmat(center(:),1,size(points,2))).^2) <= radius^2 );
        end
    end %public methods
    methods (Static=true)
        function testRand(dim, npoints, nqueries, dtype)
            %
            % test on a random sample
            % In:
            %  dim - number of dimensions
            %  npoints - number of points to sample
            %  nquerries - number of querries to test
            %  dtype - datatype to use as a handle to the casting function.
            %               (@single,@double, or @int32) are possible options
            if nargin==2
                dtype=@double;
            end
            points=dtype(rand(dim, npoints)*100);
            T=KDTree(points);
            for k=1:nqueries
                ps=sort(dtype(rand(1,dim*2)*100));
                min_corner=ps(1:dim);
                max_corner=ps(dim+1:2*dim);
                Q=T.query(min_corner, max_corner);
                KDTree.checkQuery(points, min_corner, max_corner, Q);
            end
        end
    end
    methods (Static=true, Access=protected)
        function checkQuery(points, min_corner, max_corner, result)
            % A helper for testRand
            dim=size(points,1);
            for d=1:dim
                points=points(:, points(d,:)>=min_corner(d)  & points(d,:)<=max_corner(d) );
            end
            fprintf('QSize:%i  TrueSize:%i\n', size(result,2), size(points,2));
            assert(all(size(points) == size(result)));
        end
    end
end %classdef


