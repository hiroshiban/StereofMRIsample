% ************************************************************
% This is an example of the stimulus parameter file for StereofMRIsample
% Please change the parameters below and check how these values affect
% the stimulus presentations.
%
% Created    : "2017-12-28 10:27:05 ban"
% Last Update: "2018-11-12 10:59:28 ban"
% ************************************************************

% sparam: stimulus generation parameters

%%% image size and disparity
sparam.fieldSize=[7,7];         % target stimulus size in deg
sparam.radius=[1,5];            % wedge [min,max] radius in deg
sparam.nwedges=4;               % the number of wedges
sparam.wedgeangle=80;           % wedge angle in deg
sparam.disparity=[12,  8,  4,  2, -2, -4, -8, -12]; % binocular disparities to be used, numel(sparam.disparity)=number_of_conditions
sparam.jitters=[-0.5,0,0.5];    % disparity jutters used at each stimulus presentation

%%% RDS parameters
sparam.dotRadius=[0.05,0.05];   % radius of RDS's white/black ovals
sparam.dotDens=2;               % deinsity of dot in RDS image (1-100)
sparam.colors=[255,0,128];      % RDS colors [dot1,dot2,background](0-255)
sparam.oversampling_ratio=2;    % oversampling_ratio for fine scale RDS images, [val]

% whether the random dots are filled in the background (zero-disparity) regions or not.
% if 1, the zero-disparity regions are masked.
sparam.skipzero_flg=0;

% whether avoiding dot overlaps and density biases as much as possible in generating RDSs.
% if 1, overlaps etc are taken into account. however, it is more time consuming.
sparam.avoid_bias_flg=0;

% whether using a mex function in generating RDSs.
% if 1, a mex function used. Please note that since the function put the most priority
% to processing speed, the generated RDS quality may not be always the best.
sparam.use_mex_flg=1;

%%% fixation period in sec before/after presenting the target stimuli
sparam.initial_fixation_time=[16,16]; %[2,2];

%%% stimulius presentation durations in sec
%
% The relationships of each duration parameters are as below.
%
%  ==========>>>==========>>>==========>>>==========>>>==========>>>==========>>> time
%              block                              block
%  |_______________________________|_______________________________|___ . . . 
%   stimulation_duration   blank    stimulation_duration   blank
%  |___________________|___________|___________________|__________| . . .
%  trial                           trial
%  |___|___|___|___|___|           |___|___|___|___|___| . . .
%  on off                          on off
%  |_|_|_|_| . . .                 |_|_|_|_| . . .

sparam.block_duration=32;       % block length (stimulation_duration + blank_duration)
sparam.stimulation_duration=16; % stimulation duration in one sparam.block_duration
sparam.trialDuration=2;         % trial length in sparam.stimulation_duration
sparam.stim_on_duration=1;      % stim on duration in each trial

sparam.blank_duration=sparam.block_duration-sparam.stimulation_duration;       % blank period after the stimulus presentation
sparam.stim_off_duration=sparam.trialDuration-sparam.stim_on_duration;         % stim off duration
sparam.trialsPerBlock=round(sparam.stimulation_duration/sparam.trialDuration); % trials per block

sparam.numRepeats=3;            % number of repetitions of a condition in one run

%%% fixation size and color
sparam.fixsize=24;              % the whole size (a circular hole) of the fixation cross in pixel
sparam.fixlinesize=[12,2];      % [height,width] of the fixation line in pixel
sparam.fixcolor=[255,255,255];

%%% background color
sparam.bgcolor=[128,128,128];

%%% RGB for background patches, [1x3] matrices
sparam.patch_size=[30,30];      % background patch size, [height,width] in pixels
sparam.patch_num=[20,20];       % the number of background patches along vertical and horizontal axis
sparam.patch_color1=[255,255,255];
sparam.patch_color2=[0,0,0];

%%% vernier task positions in pixels
sparam.verniersize=[6,2];       % [height,width] of the vernier bar in pixel
sparam.vernierpos=[-3,-2,-1,0,1,2,3]; % vernier position against the center of the screen in pixel

%%% loading size parameters
run(fullfile(fileparts(mfilename('fullpath')),'size_params'));
