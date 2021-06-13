**************************************************
README.txt on StereofMRIsample

Created    : "2017-12-29 13:08:04 ban"
Last Update: "2021-06-13 21:52:06 ban"
**************************************************


A block-design fMRI experiment with wedge-shaped near/far depth stimuli consisted of Rndom-Dot-Stereograms (RDSs).
function StereofMRIsample(subjID,acq,:displayfile,:stimulusfile,:gamma_table,:overwrite_flg,:force_proceed_flag,:stim_mode)
(: is optional)

[about]
- This is a sample MATLAB Psychtoolbox-3 (PTB3) script for a block-design fMRI experiment on 3D vision.
- Near/far (related to the central fixation plane) wedge stimuli defined by binocular disparities are presented.
- An example of protocol (you can change the parameters by modifying the displayfile and stimulusfile).
  16s fixation + {(16s stimulus presentation (+ 16s fixation) ) x N depth types x M cycles} + 16s fixation.
- Task during fMRIs: vernier left/right discriminations on the central fixation point.
- If you set stim_mode=2, you can present random depth patches jsut for fun.

[NOTEs]
- This script shoud be run with PTB3 or above. The script is incompatible with PTB2.
- The stimuli are created in this script in real-time with MATLAB functions in ../Generation & ../Common directories.
- displayfile & stimulusfile should be prepared as ./subjects/(subjID)/*_display.m and ./subjects/(subjID)/*_stimuli.m
  for each participant in advance of running this script.
- About the task: Attention-demanding vernier bar left/right discrimination task on the central fixation during scanning sessions
  press key1 when the vernier line is presented on the left side of the fixation (defined in displayfile)
  press key2 when the vernier line is presented on the right side of the fixation
- In general, if we use PTB3, RDS stimuli can be easily generated with Screen('DrawDots') function. However, the dots
  generated with the simple PTB3 function are not antialiased, which may cause some problem due to round-offs of the
  fine depth structures. Therefore, in this function, I am taking a different strategy to generate RDSs by putting
  antialiased (Gaussian-smoothed) dots with alpha-channel (transparency) setups and by oversampling the position shift
  (horizontal binocular disparity). That is why the stimulus generation pipeline in this function is a bit complex.

[how to run the script]
1. On the MATLAB shell, please change the working directory to
   ~/StereofMRIsample/Presentation/
2. Run the "run_exp" script
   >> run_exp('subj_name',1,1); % to present wedge-shaped depth planes.
   >> run_exp('subj_name',2,1); % to present random depth patches.
   Here, the first input variable is subject name or ID, such as 'HB' or 's01',
   the second variable should be 1 or 2,
   the third variable is run number, 1,2,3,...

For more details, please see the documents in StereofMRIsample.m
Also please see the parameter files in ~/StereofMRIsample/Presentation/subj/TEST.

[input]
sujID         : ID of a subject, a string, e.g. 'HB' or 's01'
                you have to create a directory ./subjects/(subj) and locate displayfile and stimulusfile there.
                !!!!!!!!!!!!!!!!!! IMPORTANT NOTE !!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!! if 'debug' (case insensitive) is included             !!!
                !!! in subjID string, this program runs as DEBUG mode;    !!!
                !!! stimulus images are saved as *.png format at ./images !!!
                !!!!!!!!!!!!!!!!!! IMPORTANT NOTE !!!!!!!!!!!!!!!!!!!!!!!!!!!
acq           : acquisition number (design file number),
                a integer, such as 1, 2, 3, ...
displayfile   : (optional) a display parameter file (*.m). e.g. 'nf_display_fmri.m'
                the file should be located in ./subjects/(subj)/
stimulusfile  : (optional) a stimulus generation parameter file (*.m). e.g. 'nf_stimulus_exp1.m'
                the file should be located in ./subjects/(subj)/
gamma_table   : (optional) table(s) of gamma-corrected video input values (Color LookupTable).
                256(8-bits) x 3(RGB) x 1(or 2,3,... when using multiple displays) matrix
                or a *.mat file specified with a relative path format. e.g. '/gamma_table/gamma1.mat'
                The *.mat should include a variable named "gamma_table" consists of a 256x3xN matrix.
                if you use multiple (more than 1) displays and set a 256x3x1 gamma-table, the same
                table will be applied to all displays. if the number of displays and gamma tables
                are different (e.g. you have 3 displays and 256x3x!2! gamma-tables), the last
                gamma_table will be applied to the second and third displays.
                if empty, normalized gamma table (repmat(linspace(0.0,1.0,256),3,1)) will be applied.
overwrite_flg : (optional) whether overwriting pre-existing result file. if 1, the previous results
                file with the same acquisition number will be overwritten by the previous one.
                if 0, the existing file will be backed-up by adding a prefix '_old' at the tail
                of the file. 0 by default.
force_proceed_flag : (optional) whether proceeding stimulus presentatin without waiting for
                the experimenter response (e.g. presesing the ENTER key) or a trigger.
                if 1, the stimulus presentation will be automatically carried on.
stim_mode     : (optional) stimulus mode. 1 or 2. this script presents wedge-shaped near/far depth
                stimuli when stim_mode=1, while it presents multiple random-depth patches when
                stim_mode=2 as another example. 1 by default.

[output variables]
no output matlab variable.

[output files]
1. an event log and behavior task result file
   stored ./subjects/(subjID)/results/(today)
   as ./subjects/(subjID)/results/(today)/(subjID)_(file_name_of_this_function)_run_(run_num).mat
2. a stimulus presentation log file
   stored ./subjects/(subjID)/results/(today)
   as ./subjects/(subjID)/results/(today)/(subjID)_(file_name_of_this_function)_run_(run_num).log
3. a design file
   stored ./subjects/(subjID)/design/(today)
   as ./subjects/(subjID)/design/(today)/(subjID)_design_run_(run_num).txt

[example]
>> StereofMRIsample('HB',1,'nf_display.m','nf_stimulus_exp1.m');

[About displayfile]
The contents of the displayfile is as below.

(an example of the displayfile)

% ************************************************************
% This is an example of the display parameter file for StereofMRIsample
% Please change the parameters below and check how these values affect
% the stimulus presentations.
%
% Created    : "2017-12-28 10:27:05 ban"
% Last Update: "2021-06-10 01:26:18 ban"
% ************************************************************

% dparam: display parameters

% display mode, one of "mono", "dual", "dualcross", "dualparallel", "cross", "parallel", "redgreen", "greenred",
% "redblue", "bluered", "shutter", "topbottom", "bottomtop", "interleavedline", "interleavedcolumn", "propixxmono", "propixxstereo"
dparam.ExpMode='cross';

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
dparam.fullscr='false';

dparam.ScrHeight=1024;
dparam.ScrWidth=1280;

%% shift the screen center position along y-axis (to prevent the occlusion of the stimuli due to the coil)
dparam.yshift=0;%30;

% whther skipping the PTB's vertical-sync signal test. if 1, the sync test is skipped
dparam.skip_sync_test=0;


[About stimulusfile]
The contents of the stimulusfile is as below.

(an example of the stimulusfile)

% ************************************************************
% This is an example of the stimulus parameter file for StereofMRIsample
% Please change the parameters below and check how these values affect
% the stimulus presentations.
%
% Created    : "2017-12-28 10:27:05 ban"
% Last Update: "2021-06-10 01:26:18 ban"
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

%%% for creating disparity & shadow
sparam.ipd=6.4;
sparam.pix_per_cm=57.1429;
sparam.vdist=65;
