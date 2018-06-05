%
% Mark J. Olah (mjo@cs.unm.edu)
%
% 8/2014
%
% Asyncronus ics file writing 
%
% We assume a flat data directory structure.  i.e., all files will be written
% to a single directory.
%
% Once you pass data to writeFile is is safe to delete it assuming that we are
% able to actually write.  We keep a local copy and delete than upon successful write.
%
% TODO: Async Error feedback/Unwritten data recovery

classdef ICSFileWriter < ICSAsyncBase
    methods
        function obj=ICSFileWriter(directory, Ncache sequence_field_width)
            obj=obj@ICSAsyncBase(directory, Ncache,sequence_field_width)
        end

        function status=scheduleWrite(obj, data,filename,sequence)
            %
            % Informs the C++ Class to schedule a write of data.
            % 
            % [in] [optional] sequence - integer
            if nargin==3
                sequence=[]
            end
            realFilename=obj.getRealFilename(filename, sequence)
            status=obj.call('scheduleWrite', data,
        end
        
%          function status=status(filename, sequence)
%          end
        
%          function retrieveData(filename,sequence)
%          end
    end

    
end
