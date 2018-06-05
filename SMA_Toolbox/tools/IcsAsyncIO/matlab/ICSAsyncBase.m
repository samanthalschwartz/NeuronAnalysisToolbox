%
% Mark J. Olah (mjo@cs.unm.edu)
%
% 8/2014
%
% Asyncronus ics file reading 
%
% We assume a flat data directory structure.  i.e., all files will be written
% to a single directory.
%
% Once you pass data to writeFile is is safe to delete it assuming that we are
% able to actually write.  We keep a local copy and delete than upon successful write.
%
% TODO: Async Error feedback/Unwritten data recovery

classdef ICSAsyncBase < ICSAsyncBase
    properties (Constant=true)
        statusMsg=['Success', 'Queued', 'ReadError', 'WriteError']
    end
    properties
        directory;
        Ncache=20;
        sequence=0;
    end

    methods
        function obj=ICSFileReader(directory, Ncache,sequence_field_width)
            obj.directory=directory;
        end
        
        function status=scheduleRead(filename, sequence)
        
        end

        function [data,status]=retrieveRead(filename, sequence)
            
        end
        
        function [data,status]=retrieveReadNoWait(filename, sequence)
            
        end
    end

    methods (Access=protected)
        function seq=getRealFilename(obj, filename, sequence)
            [~,f,ext]=fileparts(filename)
            if isempty(ext) || strcmp(.ics
            if isempty(sequence)
                
        end
        function d
    end
end
