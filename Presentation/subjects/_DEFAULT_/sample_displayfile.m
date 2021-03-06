% ************************************************************
% This is an example of the display parameter file for StereofMRIsample
% Please change the parameters below and check how these values affect
% the stimulus presentations.
%
% Created    : "2017-12-28 10:27:05 ban"
% Last Update: "2018-11-22 20:18:07 ban"
% ************************************************************

% dparam: display parameters

% display mode, one of "mono", "dual", "dualparallel", "dualcross", "cross", "parallel", "redgreen", "greenred",
% "redblue", "bluered", "shutter", "topbottom", "bottomtop", "interleavedline", "interleavedcolumn"
dparam.ExpMode='cross';%'dualparallel'

dparam.scrID=1; % screen ID, generally 0 for a single display setup, 1 for dual display setup

% a method to start stimulus presentation
% 0:ENTER/SPACE, 1:Left-mouse button, 2:the first MR trigger pulse (CiNet),
% 3:waiting for a MR trigger pulse (BUIC) -- checking onset of pin #11 of the parallel port,
% or 4:custom key trigger (wait for a key input that you specify as tgt_key).
dparam.start_method=4;

% a pseudo trigger key from the MR scanner when it starts, only valid when dparam.start_method=4;
dparam.custom_trigger=KbName(84); % 't' is a default trigger code from MR scanner at CiNet

%% keyboard settings
dparam.Key1=82; % key 1 'r'
dparam.Key2=71; % key 2 'g'

% otherwise, set default variables
dparam.fullscr=false;

dparam.ScrHeight=1024;
dparam.ScrWidth=1280;

%% shift the screen center position along y-axis (to prevent the occlusion of the stimuli due to the coil)
dparam.yshift=0;%30;

% whther skipping the PTB's vertical-sync signal test. if 1, the sync test is skipped
dparam.skip_sync_test=0;
