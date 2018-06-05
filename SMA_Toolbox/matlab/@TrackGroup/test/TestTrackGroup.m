
classdef TestTrackGroup < matlab.unittest.TestCase
    properties
        data_dir;
        spthsi_file;
        spthsi;
        tgrp_orig;
        tgrp;
    end

    methods (TestClassSetup)
        function loadSPTHSI(tc)
            %Get test data file
            tc.data_dir=getenv('HSM_TEST_DATA');
            if strcmp(tc.data_dir,'')
                stk = dbstack('-completenames');
                test_dir=fileparts(stk(end).file);
                tc.data_dir=fullfile(test_dir,'TestData');
            end
            if ~isdir(tc.data_dir)
                error('TestSMInteraction:construct','Test Data directroy does not exist %s', tc.data_dir);
            end
            tc.spthsi_file=fullfile(tc.data_dir, 'test.spthsi');
            
            %Load SPTHSI object 
            warning('OFF', 'SPTloadFile:dataFileNotLinked')
            tc.spthsi=SPTHSI(tc.spthsi_file);

            %Use SPTHSI object to make a 
            tc.tgrp_orig=TrackGroup(tc.spthsi);
         end
    end
    
    methods(TestMethodSetup)
        function copyTrackGroup(tc)
            tc.tgrp=tc.tgrp_orig.copy;
        end
    end

    methods (Test)
        function testConstuctor(tc)
            tc.verifyClass(tc.tgrp,'TrackGroup');
            tc.verifyNumElements(tc.tgrp.T, numel(tc.spthsi.Tracks),'Tracks copied over');
            tc.verifySize(tc.tgrp.Localizations, [sum(tc.tgrp.trackLengths()) 10],'L size');
            tc.verifyEqual(tc.tgrp.nLocalizations, size(tc.tgrp.Localizations,1),'nL');
            lengths=tc.tgrp.trackLengths();
            tc.verifyEqual(lengths,sort(lengths,'descend'), 'sorted'); 
        end
        
        %% Test Track Staistics Methods %%
        function test_trackLengths(tc)
            lengths=tc.tgrp.trackLengths();
            tc.verifySize(lengths,[tc.tgrp.nT 1]);
            tc.verifyGreaterThanOrEqual(lengths,0, 'positive');
            tc.verifyEqual(sum(lengths),tc.tgrp.nLocalizations,'sum correct');
        end

        %Also tests groupTimeSpans
        function test_trackTimeSpans(tc) 
            ts=tc.tgrp.trackTimeSpans();
            gts=tc.tgrp.groupTimeSpan();
            tc.verifySize(ts,[tc.tgrp.nT 2]);
            tc.verifyGreaterThanOrEqual(ts(:,1),0, 'start>=0');
            tc.verifyLessThanOrEqual(ts(:,1),ts(:,2), 'start<=stop');
            tc.verifyGreaterThanOrEqual(ts(:,1),gts(1), 'start>=group_start');
            tc.verifyLessThanOrEqual(ts(:,2),gts(2), 'stop<=group_stop');
        end
        
        %Also tests groupDuration
        function testDurations(tc)
            dur=tc.tgrp.trackDurations();
            gdur=tc.tgrp.groupDuration();
            gts=tc.tgrp.groupTimeSpan();
            ts=tc.tgrp.trackTimeSpans();
            tc.verifySize(dur,[tc.tgrp.nT 1], 'size');
            tc.verifyGreaterThanOrEqual(dur(:),0, '>=0');
            tc.verifyEqual(dur, ts(:,2)-ts(:,1), 'durations match times');
            tc.verifyGreaterThanOrEqual(gdur, max(dur), 'group_duration>=track_duration');
            tc.verifyEqual(gdur, gts(2)-gts(1), 'group_duration matched group_times');
        end

        %Also tests groupSpectraBounds
        function test_trackSpectraBounds(tc)
            bounds=tc.tgrp.trackSpectraBounds();
            gbounds=tc.tgrp.groupSpectraBound();
            tc.verifySize(bounds,[tc.tgrp.nT 2]);
            tc.verifySize(gbounds,[1 2]);
            tc.verifyGreaterThanOrEqual(bounds,100, 'lambda>100');
            tc.verifyLessThanOrEqual(bounds,1000, 'lambda<1000');
            tc.verifyGreaterThanOrEqual(gbounds,100, 'group_lambda>100');
            tc.verifyLessThanOrEqual(gbounds,1000, 'group_lambda<1000');
            tc.verifyLessThanOrEqual(bounds(:,1),bounds(:,2), 'ordered');
            tc.verifyLessThanOrEqual(gbounds(1),gbounds(2), 'group ordered');
            tc.verifyGreaterThanOrEqual(bounds(:,1),gbounds(1), 'group lower bound');
            tc.verifyLessThanOrEqual(gbounds(:,2),gbounds(2), 'group upper bound');
        end

        %Also tests groupBBox
        function test_trackBBoxes(tc)
            bboxes=tc.tgrp.trackBBoxes();
            gbbox=tc.tgrp.groupBBox();
            tc.verifySize(bboxes,[tc.tgrp.nT 4]);
            tc.verifyLessThanOrEqual(bboxes(:,1),bboxes(:,3), 'xmin<=xmax');
            tc.verifyLessThanOrEqual(bboxes(:,2),bboxes(:,4), 'ymin<=ymax');
            tc.verifyGreaterThanOrEqual(bboxes(:,1),gbbox(:,1), 'xmin>=group_xmin');
            tc.verifyGreaterThanOrEqual(bboxes(:,2),gbbox(:,2), 'ymin>=group_ymin');
            tc.verifyLessThanOrEqual(bboxes(:,3),gbbox(:,3), 'xmax<=group_xmax');
            tc.verifyLessThanOrEqual(bboxes(:,4),gbbox(:,4), 'ymax<=group_ymax');
        end
       
        function testFilterLength(tc)
            null_filter_val=0;
            filter_val=10;
            mytgrp=tc.tgrp;
            tc.testFilter(@mytgrp.filterLength, null_filter_val, filter_val);
            lengths=tc.tgrp.trackLengths();
            tc.verifyGreaterThanOrEqual(lengths,filter_val,'filtered tracks removed');
        end

        function testFilterDuration(tc)
            null_filter_val=0;
            filter_val=0.1;
            mytgrp=tc.tgrp;
            tc.testFilter(@mytgrp.filterDuration, null_filter_val, filter_val);
            dur=tc.tgrp.trackDurations();
            tc.verifyGreaterThanOrEqual(dur,filter_val,'filtered tracks removed');
        end

        function testFilterSpectralDeviation(tc)
            gbound=tc.tgrp.groupSpectraBound();
            null_filter_val=gbound(2)-gbound(1);
            filter_val=15;
            mytgrp=tc.tgrp;
            tc.testFilter(@mytgrp.filterSpectralDeviation, null_filter_val, filter_val);
            bounds=tc.tgrp.trackSpectraBounds();
            tc.verifyLessThanOrEqual(bounds(:,2)-bounds(:,1),filter_val,'filtered tracks removed');
        end

        function testFilterSelection(tc)
            null_filter_val=[];
            filter_val=1:10;
            mytgrp=tc.tgrp;
            tc.testFilter(@mytgrp.filterSelection, null_filter_val, filter_val);
        end
        
        function testFilterRandomSubsample(tc)
            null_filter_val=1.0;
            filter_val=0.25;
            old_nT=tc.tgrp.nT;
            mytgrp=tc.tgrp;
            tc.testFilter(@mytgrp.filterRandomSubsample, null_filter_val, filter_val);
            tc.verifyEqual(tc.tgrp.nT,round(filter_val*old_nT),'correct fraction removed');
        end

    end

    methods
        function testFilter(tc, filter_func, null_filter_val, filter_val)
            %Test that null_filter_val does nothing
            old_nT=tc.tgrp.nT;
            nrejected=filter_func(null_filter_val);
            tc.verifyEqual(tc.tgrp.nT,old_nT,'null val leaves tracks unchanged');
            tc.verifyEqual(nrejected,0,'nrejected 0');

            %Test that filter_val does something
            old_nT=tc.tgrp.nT;
            nrejected=filter_func(filter_val);
            tc.verifyLessThan(tc.tgrp.nT,old_nT,'filter removes tracks');
            tc.verifyEqual(nrejected,old_nT-tc.tgrp.nT,'nrejected accurate');
        end
    end 
end
    
