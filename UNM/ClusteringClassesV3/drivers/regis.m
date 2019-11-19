% Perform registration analysis use Mark Olah's RegistrationAnalysis class.
% Required: quadView_SplitDIPimage.m (or quadView_SplitImage.m)
% Michael Wester (12/5/16)

% Mark Olah's RegistrationAnalysis software needs to be on the MATLAB path:
%addpath(genpath('~/Documents/MATLAB/SMA_Toolbox/development/matlab'));

% The top right and bottom left are the channels that correspond to
% A488 and A647 data, respectively.  A background image is also needed.
% The user needs to provide:
%    Pixel Size [X, y] (um)
%    CCD Gain (e-/ADU)
%    CCD Background (ADU)
% The latter two (CCD Gain and Background) can also be computed in the GUI via
% Calibrate.

% Split up the 4 channels as DIPimage images.
iodir = 'mat';
load(fullfile(iodir, 'GRID.mat'));
[BL, TL, TR, BR] = quadView_SplitDIPimage(Data);
%[BL, TL, TR, BR] = quadView_SplitImage(Data);

% Assign the quadrants as appropriate.
A488 = TR(:, :, 50:99);
A647 = BL(:, :, 50:99);
BG   = BR(:, :, 50:99); % fake for now

% Save the sequences in files for further processing by RegistrationAnalysis.
% The first three files are only used by Calibrate.
sequence = A488;
save(fullfile(iodir, 'Beads488'), 'sequence');
sequence = A647;
save(fullfile(iodir, 'Beads647'), 'sequence');
sequence = BG;
save(fullfile(iodir, 'Background'), 'sequence');
% A488_647 has the 488 channel on the left, the 647 channel on the right.
% This file is used to determine the registration.
sequence = [A488, A647];
save(fullfile(iodir, 'A488_647'), 'sequence');
% A647_488 has the 647 channel on the left, the 488 channel on the right.
%sequence = [A647, A488];
%save(fullfile(iodir, 'A647_488'), 'sequence');

% Some typical values:
%    Pixel Size [X, Y](um): [0.104, 0.104]
%    CCD Gain (e-/ADU):     (488) 0.57275  (647) 0.059015
%    CCD Background (ADU):  (488,647) 102.5624
%
% Important Advanced Params:
%    Boundary Width
%    Min Intensity
%    Max Reference Displacement (2-4 pixels)
%
% Notes:
%    .reganalysis file is self-contained
%    Reprocess ALL Sequences => Reprocess ALL Maps
%    RA.data, RA.maps contain useful results

RA = RegistrationAnalysis();

RA.pixelSize = [0.104, 0.104];
RA.CCDBackground = 102.5264;
RA.CCDGain = 0.059015;
RA.workingDir = iodir;
RA.Params.MinIntensity = 200;

RA.load(fullfile(iodir, 'A488_647.mat'));
%RA.load(fullfile(iodir, 'A647_488.mat'));
% Start the GUI.
RA.gui();

% Save the RegistrationAnalysis object ...
% ... in Mark's .reganalysis file,
RA.save();
% ... as a regular MATLAB save file (for use by run_RA_PC_SR.m).
save(fullfile(iodir, 'RegAnal.mat'), 'RA');

% The mapping takes the right channel into the left.  Note that mapping occurs
% in absolute coordinates, so the right channel will have larger x-coordinates
% than the left.  This becomes important if data from the two channels is
% imported separately in _relative_ coordinates, in which case the
% x-coordinates of the isolated right channel must be adjusted to use the
% mapping properly.  See run_RA_PC_SR.m
%
% Notes from Samantha Schwartz:
% Also, from a tracked RPT object (rpt) and a saved RegistrationAnalysis
% object (RA) you can call the following to get the shifted tracks:
% tracks = RA.shiftRPTTracksPixels(rpt);
% The tracks returned will be in the same format as they are in the RPT object.
