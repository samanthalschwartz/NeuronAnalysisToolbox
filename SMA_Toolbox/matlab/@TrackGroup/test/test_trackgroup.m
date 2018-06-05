%This runs all the tests in the TestSMInteraction class
%if data_dir is given, pass this on through the environment variable
%HSM_TEST_DATA
function result=test_trackgroup(data_dir)

    import matlab.unittest.TestSuite;

    if nargin==1
        setenv('HSM_TEST_DATA',data_dir)
    end

    suite = TestSuite.fromClass(?TestTrackGroup);
    result =run(suite);

end
