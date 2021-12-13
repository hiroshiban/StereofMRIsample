function StereofMRIsample(subjID,acq,displayfile,stimulusfile,gamma_table,overwrite_flg,force_proceed_flag,stim_mode)

% A block-design fMRI experiment with wedge-shaped near/far depth stimuli consisted of Rndom-Dot-Stereograms (RDSs).
% function StereofMRIsample(subjID,acq,:displayfile,:stimulusfile,:gamma_table,:overwrite_flg,:force_proceed_flag,:stim_mode)
% (: is optional)
%
% [about]
% - This is a sample MATLAB Psychtoolbox-3 (PTB3) script for a block-design fMRI experiment on 3D vision.
% - Near/far (related to the central fixation plane) wedge stimuli defined by binocular disparities are presented.
% - An example of protocol (you can change the parameters by modifying the displayfile and stimulusfile).
%   16s fixation + {(16s stimulus presentation (+ 16s fixation) ) x N depth types x M cycles} + 16s fixation.
% - Task during fMRIs: vernier left/right discriminations on the central fixation point.
% - If you set stim_mode=2, you can present random depth patches jsut for fun.
%
% [NOTEs]
% - This script shoud be run with PTB3 or above. The script is incompatible with PTB2.
% - The stimuli are created in this script in real-time with MATLAB functions in ../Generation & ../Common directories.
% - displayfile & stimulusfile should be prepared as ./subjects/(subjID)/*_display.m and ./subjects/(subjID)/*_stimuli.m
%   for each participant in advance of running this script.
% - About the task: Attention-demanding vernier bar left/right discrimination task on the central fixation during scanning sessions
%   press key1 when the vernier line is presented on the left side of the fixation (defined in displayfile)
%   press key2 when the vernier line is presented on the right side of the fixation
% - In general, if we use PTB3, RDS stimuli can be easily generated with Screen('DrawDots') function. However, the dots
%   generated with the simple PTB3 function are not antialiased, which may cause some problem due to round-offs of the
%   fine depth structures. Therefore, in this function, I am taking a different strategy to generate RDSs by putting
%   antialiased (Gaussian-smoothed) dots with alpha-channel (transparency) setups and by oversampling the position shift
%   (horizontal binocular disparity). That is why the stimulus generation pipeline in this function is a bit complex.
%
% [how to run the script]
% 1. On the MATLAB shell, please change the working directory to
%    ~/StereofMRIsample/Presentation/
% 2. Run the "run_exp" script
%    >> run_exp('subj_name',1,1); % to present wedge-shaped depth planes.
%    >> run_exp('subj_name',2,1); % to present random depth patches.
%    Here, the first input variable is subject name or ID, such as 'HB' or 's01',
%    the second variable should be 1 or 2,
%    the third variable is run number, 1,2,3,...
%
% Created    : "2017-12-29 14:33:31 ban"
% Last Update: "2021-12-14 01:13:09 ban"
%
%
% [input]
% sujID         : ID of a subject, a string, e.g. 'HB' or 's01'
%                 you have to create a directory ./subjects/(subj) and locate displayfile and stimulusfile there.
%                 !!!!!!!!!!!!!!!!!! IMPORTANT NOTE !!!!!!!!!!!!!!!!!!!!!!!!
%                 !!! if 'debug' (case insensitive) is included          !!!
%                 !!! in subjID string, this program runs as DEBUG mode; !!!
%                 !!! stimulus images are saved as *.png format at       !!!
%                 !!! ~/StereofMRIsample/Presentation/images             !!!
%                 !!!!!!!!!!!!!!!!!! IMPORTANT NOTE !!!!!!!!!!!!!!!!!!!!!!!!
% acq           : acquisition number (design file number),
%                 a integer, such as 1, 2, 3, ...
% displayfile   : (optional) a display parameter file (*.m). e.g. 'nf_display_fmri.m'
%                 the file should be located in ./subjects/(subj)/
% stimulusfile  : (optional) a stimulus generation parameter file (*.m). e.g. 'nf_stimulus_exp1.m'
%                 the file should be located in ./subjects/(subj)/
% gamma_table   : (optional) table(s) of gamma-corrected video input values (Color LookupTable).
%                 256(8-bits) x 3(RGB) x 1(or 2,3,... when using multiple displays) matrix
%                 or a *.mat file specified with a relative path format. e.g. '/gamma_table/gamma1.mat'
%                 The *.mat should include a variable named "gamma_table" consists of a 256x3xN matrix.
%                 if you use multiple (more than 1) displays and set a 256x3x1 gamma-table, the same
%                 table will be applied to all displays. if the number of displays and gamma tables
%                 are different (e.g. you have 3 displays and 256x3x!2! gamma-tables), the last
%                 gamma_table will be applied to the second and third displays.
%                 if empty, normalized gamma table (repmat(linspace(0.0,1.0,256),3,1)) will be applied.
% overwrite_flg : (optional) whether overwriting pre-existing result file. if 1, the previous results
%                 file with the same acquisition number will be overwritten by the previous one.
%                 if 0, the existing file will be backed-up by adding a prefix '_old' at the tail
%                 of the file. 0 by default.
% force_proceed_flag : (optional) whether proceeding stimulus presentatin without waiting for
%                 the experimenter response (e.g. presesing the ENTER key) or a trigger.
%                 if 1, the stimulus presentation will be automatically carried on.
% stim_mode     : (optional) stimulus mode. 1 or 2. this script presents wedge-shaped near/far depth
%                 stimuli when stim_mode=1, while it presents multiple random-depth patches when
%                 stim_mode=2 as another example. 1 by default.
%
% [output variables]
% no output matlab variable.
%
% [output files]
% 1. an event log and behavior task result file
%    stored ./subjects/(subjID)/results/(today)
%    as ./subjects/(subjID)/results/(today)/(subjID)_(file_name_of_this_function)_run_(run_num).mat
% 2. a stimulus presentation log file
%    stored ./subjects/(subjID)/results/(today)
%    as ./subjects/(subjID)/results/(today)/(subjID)_(file_name_of_this_function)_run_(run_num).log
% 3. a design file
%    stored ./subjects/(subjID)/design/(today)
%    as ./subjects/(subjID)/design/(today)/(subjID)_design_run_(run_num).txt
%
% [example]
% >> StereofMRIsample('HB',1,'nf_display.m','nf_stimulus_exp1.m');
%
% [About displayfile]
% The contents of the displayfile is as below.
%
% (an example of the displayfile)
%
% % ************************************************************
% % This is an example of the display parameter file for StereofMRIsample
% % Please change the parameters below and check how these values affect
% % the stimulus presentations.
% %
% % Created    : "2017-12-28 10:27:05 ban"
% % Last Update: "2021-06-10 01:26:18 ban"
% % ************************************************************
%
% % dparam: display parameters
%
% % display mode, one of "mono", "dual", "dualcross", "dualparallel", "cross", "parallel", "redgreen", "greenred",
% % "redblue", "bluered", "shutter", "topbottom", "bottomtop", "interleavedline", "interleavedcolumn", "propixxmono", "propixxstereo"
% dparam.ExpMode='cross';
%
% dparam.scrID=1; % screen ID, generally 0 for a single display setup, 1 for dual display setup
%
% % a method to start stimulus presentation
% % 0:ENTER/SPACE, 1:Left-mouse button, 2:the first MR trigger pulse (CiNet),
% % 3:waiting for a MR trigger pulse (BUIC) -- checking onset of pin #11 of the parallel port,
% % or 4:custom key trigger (wait for a key input that you specify as tgt_key).
% dparam.start_method=4;
%
% % a pseudo trigger key from the MR scanner when it starts, only valid when dparam.start_method=4;
% dparam.custom_trigger=KbName(84); % 't' is a default trigger code from MR scanner at CiNet
%
% %% keyboard settings
% dparam.Key1=82; % key 1 'r'
% dparam.Key2=71; % key 2 'g'
%
% % otherwise, set default variables
% dparam.fullscr='false';
%
% dparam.ScrHeight=1024;
% dparam.ScrWidth=1280;
%
% %% shift the screen center position along y-axis (to prevent the occlusion of the stimuli due to the coil)
% dparam.yshift=0;%30;
%
% % whther skipping the PTB's vertical-sync signal test. if 1, the sync test is skipped
% dparam.skip_sync_test=0;
%
%
% [About stimulusfile]
% The contents of the stimulusfile is as below.
%
% (an example of the stimulusfile)
%
% % ************************************************************
% % This is an example of the stimulus parameter file for StereofMRIsample
% % Please change the parameters below and check how these values affect
% % the stimulus presentations.
% %
% % Created    : "2017-12-28 10:27:05 ban"
% % Last Update: "2021-06-10 01:26:18 ban"
% % ************************************************************
%
% % sparam: stimulus generation parameters
%
% %%% image size and disparity
% sparam.fieldSize=[7,7];         % target stimulus size in deg
% sparam.radius=[1,5];            % wedge [min,max] radius in deg
% sparam.nwedges=4;               % the number of wedges
% sparam.wedgeangle=80;           % wedge angle in deg
% sparam.disparity=[12,  8,  4,  2, -2, -4, -8, -12]; % binocular disparities to be used, numel(sparam.disparity)=number_of_conditions
% sparam.jitters=[-0.5,0,0.5];    % disparity jutters used at each stimulus presentation
%
% %%% RDS parameters
% sparam.dotRadius=[0.05,0.05];   % radius of RDS's white/black ovals
% sparam.dotDens=2;               % deinsity of dot in RDS image (1-100)
% sparam.colors=[255,0,128];      % RDS colors [dot1,dot2,background](0-255)
% sparam.oversampling_ratio=2;    % oversampling_ratio for fine scale RDS images, [val]
%
% % whether the random dots are filled in the background (zero-disparity) regions or not.
% % if 1, the zero-disparity regions are masked.
% sparam.skipzero_flg=0;
%
% % whether avoiding dot overlaps and density biases as much as possible in generating RDSs.
% % if 1, overlaps etc are taken into account. however, it is more time consuming.
% sparam.avoid_bias_flg=0;
%
% % whether using a mex function in generating RDSs.
% % if 1, a mex function used. Please note that since the function put the most priority
% % to processing speed, the generated RDS quality may not be always the best.
% sparam.use_mex_flg=1;
%
% %%% fixation period in sec before/after presenting the target stimuli
% sparam.initial_fixation_time=[16,16]; %[2,2];
%
% %%% stimulius presentation durations in sec
% %
% % The relationships of each duration parameters are as below.
% %
% %  ==========>>>==========>>>==========>>>==========>>>==========>>>==========>>> time
% %              block                              block
% %  |_______________________________|_______________________________|___ . . .
% %   stimulation_duration   blank    stimulation_duration   blank
% %  |___________________|___________|___________________|__________| . . .
% %  trial                           trial
% %  |___|___|___|___|___|           |___|___|___|___|___| . . .
% %  on off                          on off
% %  |_|_|_|_| . . .                 |_|_|_|_| . . .
%
% sparam.block_duration=32;       % block length (stimulation_duration + blank_duration)
% sparam.stimulation_duration=16; % stimulation duration in one sparam.block_duration
% sparam.trialDuration=2;         % trial length in sparam.stimulation_duration
% sparam.stim_on_duration=1;      % stim on duration in each trial
%
% sparam.blank_duration=sparam.block_duration-sparam.stimulation_duration;       % blank period after the stimulus presentation
% sparam.stim_off_duration=sparam.trialDuration-sparam.stim_on_duration;         % stim off duration
% sparam.trialsPerBlock=round(sparam.stimulation_duration/sparam.trialDuration); % trials per block
%
% sparam.numRepeats=3;            % number of repetitions of a condition in one run
%
% %%% fixation size and color
% sparam.fixsize=24;              % the whole size (a circular hole) of the fixation cross in pixel
% sparam.fixlinesize=[12,2];      % [height,width] of the fixation line in pixel
% sparam.fixcolor=[255,255,255];
%
% %%% background color
% sparam.bgcolor=[128,128,128];
%
% %%% RGB for background patches, [1x3] matrices
% sparam.patch_size=[30,30];      % background patch size, [height,width] in pixels
% sparam.patch_num=[20,20];       % the number of background patches along vertical and horizontal axis
% sparam.patch_color1=[255,255,255];
% sparam.patch_color2=[0,0,0];
%
% %%% vernier task positions in pixels
% sparam.verniersize=[6,2];       % [height,width] of the vernier bar in pixel
% sparam.vernierpos=[-3,-2,-1,0,1,2,3]; % vernier position against the center of the screen in pixel
%
% %%% for creating disparity & shadow
% sparam.ipd=6.4;
% sparam.pix_per_cm=57.1429;
% sparam.vdist=65;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Check the input variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear global; clear mex;
if nargin<2, help(mfilename()); return; end
if nargin<3 || isempty(displayfile), displayfile=[]; end
if nargin<4 || isempty(stimulusfile), stimulusfile=[]; end
if nargin<5 || isempty(gamma_table), gamma_table=[]; end
if nargin<6 || isempty(overwrite_flg), overwrite_flg=0; end
if nargin<7 || isempty(force_proceed_flag), force_proceed_flag=0; end
if nargin<8 || isempty(stim_mode), stim_mode=1; end

% check the aqcuisition number.
if acq<1, error('Acquistion number must be integer and greater than zero'); end

rootDir=fileparts(mfilename('fullpath'));

% check the subject directory
if ~exist(fullfile(pwd,'subjects',subjID),'dir'), error('can not find subj directory. check the input variable.'); end

% check the display/stimulus files
if ~isempty(displayfile)
  if ~strcmpi(displayfile(end-1:end),'.m'), displayfile=[displayfile,'.m']; end
  if ~exist(fullfile(rootDir,'subjects',subjID,displayfile),'file'), error('displayfile not found. check the input variable.'); end
end

if ~isempty(stimulusfile)
  if ~strcmpi(stimulusfile(end-1:end),'.m'), stimulusfile=[stimulusfile,'.m']; end
  if ~exist(fullfile(rootDir,'subjects',subjID,stimulusfile),'file'), error('stimulusfile not found. check the input variable.'); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Add paths to the subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add paths to the subfunctions
addpath(genpath(fullfile(rootDir,'..','Common')));
addpath(fullfile(rootDir,'..','Generation'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% For a log file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% result directry & file
resultDir=fullfile(rootDir,'subjects',num2str(subjID),'results',datestr(now,'yymmdd'));
if ~exist(resultDir,'dir'), mkdir(resultDir); end

% record the output window
logfname=fullfile(resultDir,[num2str(subjID),sprintf('_%s_run_',mfilename()),num2str(acq,'%02d'),'.log']);
diary(logfname);
warning off; %#ok warning('off','MATLAB:dispatcher:InexactCaseMatch');


%%%%% try & catch %%%%%
try


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Check the PTB version
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PTB_OK=CheckPTBversion(3); % check wether the PTB version is 3
if ~PTB_OK, error('Wrong version of Psychtoolbox is running. %s requires PTB ver.3',mfilename()); end

% debug level, black screen during calibration
Screen('Preference','VisualDebuglevel',3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Setup random seed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

InitializeRandomSeed();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Reset display Gamma-function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(gamma_table)
  gamma_table=repmat(linspace(0.0,1.0,256),3,1)'; %#ok
  GammaResetPTB(1.0);
else
  GammaLoadPTB(gamma_table);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Validate dparam (displayfile) and sparam (stimulusfile) structures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% organize dparam
dparam=struct(); % initialize
if ~isempty(displayfile), run(fullfile(rootDir,'subjects',subjID,displayfile)); end % load specific dparam parameters configured for each of the participants
dparam=ValidateStructureFields(dparam,... % validate fields and set the default values to missing field(s)
         'ExpMode','cross',...
         'scrID',0,...
         'start_method',1,...
         'custom_trigger',KbName(84),...
         'Key1',37,...
         'Key2',39,...
         'fullscr',false,...
         'ScrHeight',1200,...
         'ScrWidth',1920,...
         'yshift',0,...
         'skip_sync_test',0);

% organize sparam
sparam=struct(); % initialize
if ~isempty(stimulusfile), run(fullfile(rootDir,'subjects',subjID,stimulusfile)); end % load specific sparam parameters configured for each of the participants
sparam=ValidateStructureFields(sparam,... % validate fields and set the default values to missing field(s)
         'fieldSize',[7,7],...
         'radius',[1,5],...
         'nwedges',4,...
         'wedgeangle',80,...
         'disparity',[12,  8,  4,  2, -2, -4, -8, -12],...
         'jitters',[-0.5,0,0.5],...
         'dotRadius',[0.05,0.05],...
         'dotDens',2,...
         'colors',[255,0,128],...
         'oversampling_ratio',2,...
         'skipzero_flg',0,...
         'avoid_bias_flg',0,...
         'use_mex_flg',1,...
         'initial_fixation_time',[16,16],...
         'block_duration',32,...
         'stimulation_duration',16,...
         'trialDuration',2,...
         'stim_on_duration',1,...
         'numRepeats',3,...
         'fixsize',24,...
         'fixlinesize',[12,2],...
         'fixcolor',[255,255,255],...
         'bgcolor',[128,128,128],...
         'patch_size',[30,30],...
         'patch_num',[20,20],...
         'patch_color1',[255,255,255],...
         'patch_color2',[0,0,0],...
         'verniersize',[6,2],...
         'vernierpos',[-3,-2,-1,0,1,2,3],...
         'ipd',6.4,...
         'pix_per_cm',57.1429,...
         'vdist',65);

% set the other parameters
dparam.RunScript=mfilename();
sparam.RunScript=mfilename();

sparam.numConds=numel(sparam.disparity); % the number of conditions

% set the trials per block during the blank period
blank_trialsPerBlock=round(sparam.blank_duration/sparam.trialDuration);

% validate the stimulus file contents, you can add more checks for your safety
if min(sparam.fieldSize)<max(sparam.radius), error('sparam.radius should be smaller than sparam.fieldSize. check the stimulusfile.'); end
if mod(sparam.stimulation_duration,sparam.trialDuration), error('sparam.stimulation_duration can not be divided by sparam.trialDuration. check the stimulusfile.'); end
if mod(sparam.blank_duration,sparam.trialDuration), error('sparam.blank_duration can not be divided by sparam.trialDuration. check the stimulusfile.'); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Displaying the presentation parameters you set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nThe Presentation Parameters are as below.\n\n');
fprintf('************************************************\n');
fprintf('****** Script, Subject, Acquistion Number ******\n');
fprintf('Running Script Name    : %s\n',mfilename());
fprintf('Subject ID             : %s\n',subjID);
fprintf('Acquisition Number     : %d\n',acq);
fprintf('********* Run Type, Display Image Type *********\n');
fprintf('Display Mode           : %s\n',dparam.ExpMode);
fprintf('use Full Screen Mode   : %d\n',dparam.fullscr);
fprintf('Start Method           : %d\n',dparam.start_method);
if dparam.start_method==4
  fprintf('Custom Trigger         : %s\n',dparam.custom_trigger);
end
fprintf('*************** Screen Settings ****************\n');
fprintf('Screen Height          : %d\n',dparam.ScrHeight);
fprintf('Screen Width           : %d\n',dparam.ScrWidth);
fprintf('*********** Stimulus Conditions etc. ***********\n');
fprintf('number of conditions   : %d\n',sparam.numConds);
fprintf('*********** Stimulation Periods etc. ***********\n');
fprintf('Fixation Time(sec)     : [%d,%d]\n',...
     sparam.initial_fixation_time(1),sparam.initial_fixation_time(2));
fprintf('Block Duration(sec)    : %d\n',sparam.block_duration);
fprintf('Trial Duration(sec)    : %d\n',sparam.stimulation_duration);
fprintf('Blank Duration(sec)    : %d\n',sparam.blank_duration);
fprintf('Repetitions(cycles)    : %d\n',sparam.numRepeats);
fprintf('Total Time (sec)       : %d\n',...
     sum(sparam.initial_fixation_time)+sparam.numConds*sparam.numRepeats*(sparam.block_duration));
fprintf('************ Response key settings *************\n');
fprintf('Reponse Key #1         : %d=%s\n',dparam.Key1,KbName(dparam.Key1));
fprintf('Reponse Key #2         : %d=%s\n',dparam.Key2,KbName(dparam.Key2));
fprintf('************************************************\n\n');
fprintf('Please carefully check before proceeding.\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Creating stimulus presentation protocol, a design structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create design array (stimulus presentation sequences)
design=GenerateRandomDesignSequence(sparam.numConds,sparam.numRepeats,1,0,1)';

% save the current design file
designDir=fullfile(rootDir,'subjects',num2str(subjID),'design',datestr(now,'yymmdd'));
if ~exist(designDir,'dir'), mkdir(designDir); end

% open the design file
fid=fopen(fullfile(designDir,[num2str(subjID),'_design_run_',num2str(acq,'%02d'),'.txt']),'w');
if fid==-1, error('can not write the design file. pleaes run fclose all and re-run the script.'); end

% adjust the design contents with padding 0 (fixation or blank period)
designcur=design;

if sparam.blank_duration~=0 % insert zero (=blank condition) between stimulation condition
  designcur=([design,zeros(numel(design),1)])';
  designcur=designcur(:);
end

if sparam.initial_fixation_time(1)~=0 % insert zero at the head of the design (initial fixation)
  if designcur(1)~=0, designcur=[0;designcur]; end
end

if sparam.initial_fixation_time(2)~=0 % insert zero at the head of the design (final fixation)
  if designcur(end)~=0, designcur=[designcur;0]; end
end

% write the stimulus presentation order (condition order) to the design file
for ii=1:1:numel(designcur), fprintf(fid,'%d\n',designcur(ii)); end
fclose(fid);

clear designcur;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Initialize response & event logger objects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize MATLAB objects for event and response logs
event=eventlogger();
resps=responselogger([dparam.Key1,dparam.Key2]);
resps.initialize(event); % initialize responselogger


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Wait for user reponse to start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~force_proceed_flag
  [user_answer,resps]=resps.wait_to_proceed();
  if ~user_answer, diary off; return; end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Initialization of Left & Right screens for binocular presenting/viewing mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dparam.skip_sync_test, Screen('Preference','SkipSyncTests',1); end

% ************************************* IMPORTANT NOTE *****************************************
% if the console PC has been connected to two 3D displays with the expanding display setups and
% some shutter goggles (e.g. nVidia 3DVision2) are used for displaying 3D stimulus with MATLAB
% Psychtoolbox3 (PTB3), the 3D stimuli can be presented properly only when we select the first
% display (scrID=1) for stimulus presentations, while left/right images seem to be flipped if
% we select the second display (scrID=2) as the main stimulus presentation window. This may be
% a bug of PTB3. So in any case, if you run 3D vision experiments with dual monitors, it would
% be safer to always chose the first monitor for stimulus presentations. Please be careful.
% ************************************* IMPORTANT NOTE *****************************************

[winPtr,winRect,nScr,dparam.fps,dparam.ifi,initDisplay_OK]=InitializePTBDisplays(dparam.ExpMode,sparam.bgcolor,0,[],dparam.scrID);
if ~initDisplay_OK, error('Display initialization error. Please check your exp_run parameter.'); end
HideCursor();

dparam.fps=60; % set the fixed flips/sec velue just in case, as the PTB sometimes underestimates the actual vertical sync signals.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Setting the PTB runnning priority to MAX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set the priority of this script to MAX
priorityLevel=MaxPriority(winPtr,'WaitBlanking');
Priority(priorityLevel);

% conserve VRAM memory: Workaround for flawed hardware and drivers
% 32 == kPsychDontShareContextRessources: Do not share ressources between
% different onscreen windows. Usually you want PTB to share all ressources
% like offscreen windows, textures and GLSL shaders among all open onscreen
% windows. If that causes trouble for some weird reason, you can prevent
% automatic sharing with this flag.
%Screen('Preference','ConserveVRAM',32);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Setting the PTB OpenGL functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Enable OpenGL mode of Psychtoolbox: This is crucially needed for clut animation
InitializeMatlabOpenGL();

% This script calls Psychtoolbox commands available only in OpenGL-based
% versions of the Psychtoolbox. (So far, the OS X Psychtoolbox is the
% only OpenGL-base Psychtoolbox.)  The Psychtoolbox command AssertPsychOpenGL will issue
% an error message if someone tries to execute this script on a computer without
% an OpenGL Psychtoolbox
AssertOpenGL();

% set OpenGL blend functions
Screen('BlendFunction',winPtr,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Initializing MATLAB OpenGL shader API
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Not required for the current display and stimulus setups
% just call DrawTextureWithCLUT with window pointer alone
%DrawTextureWithCLUT(winPtr);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Displaying 'Initializing...'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% displaying texts on the center of the screen
%DisplayMessage2('Initializing...',sparam.bgcolor,winPtr,nScr,'Arial',36);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Initializing variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% unit conversions

% cm per pix
sparam.cm_per_pix=1/sparam.pix_per_cm;

% pixles per degree
sparam.pix_per_deg=round( 1/( 180*atan(sparam.cm_per_pix/sparam.vdist)/pi ) );

% calculate wedge height
wedge_height=zeros(numel(sparam.jitters),numel(sparam.disparity));
for tt=1:1:numel(sparam.jitters) % for sparam.jitter, e.g. -1,0,+1
  for ii=1:1:numel(sparam.disparity)
    wedge_height(tt,ii)=CalcDistFromDisparity(sparam.ipd,sparam.disparity(ii)+sparam.jitters(tt),sparam.vdist); % unit: cm;
  end
end

% generate base wedge field
wedge_field=cell(sparam.trialsPerBlock,1);
for ii=1:1:sparam.trialsPerBlock
  wedge_field{ii}=nf_CreateWedgesField(sparam.fieldSize,sparam.radius(1),sparam.radius(2),1,sparam.nwedges,sparam.wedgeangle,...
                                       -1*sparam.wedgeangle/2+(ii-1)*(90-sparam.wedgeangle/2)/sparam.trialsPerBlock,...
                                       sparam.pix_per_deg,sparam.oversampling_ratio);
end

% adjust parameters for oversampling (this adjustment shoud be done after creating heightfields)
dotDens=sparam.dotDens/sparam.oversampling_ratio;
ipd=sparam.ipd*sparam.oversampling_ratio;
vdist=sparam.vdist*sparam.oversampling_ratio;
pix_per_cm_x=sparam.pix_per_cm*sparam.oversampling_ratio;
pix_per_cm_y=sparam.pix_per_cm;

% generate ovals to be used in RDS
dotSize=round(sparam.dotRadius.*[pix_per_cm_y,pix_per_cm_x]*2); % radius(cm) --> diameter(pix)
basedot=double(MakeFineOval(dotSize,[sparam.colors,0],sparam.colors(3),1.2,2,1,0,0));
wdot=basedot(:,:,1);     % get only gray scale image (white)
bdot=basedot(:,:,2);     % get only gray scale image (black)
dotalpha=basedot(:,:,4)./max(max(basedot(:,:,4))); % get alpha channel value 0-1.0;

% calculate position shifts for each value of heightfield
posL=cell(numel(sparam.jitters),1);
posR=cell(numel(sparam.jitters),1);
for tt=1:1:numel(sparam.jitters)
  posL{tt}=cell(numel(sparam.disparity),sparam.trialsPerBlock);
  posR{tt}=cell(numel(sparam.disparity),sparam.trialsPerBlock);
end

for tt=1:1:numel(sparam.jitters) % for jitter, -1,0,+1
  for ii=1:1:numel(sparam.disparity)
    for jj=1:1:sparam.trialsPerBlock
      [posL{tt}{ii,jj},posR{tt}{ii,jj}]=RayTrace_ScreenPos_X_MEX(wedge_height(tt,ii)*wedge_field{jj},ipd,vdist,pix_per_cm_x,0);
    end
  end
end
clear wedge_field wedge_height;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Debug codes
%%%% saving the stimulus images as *.png format files and enter the debug (keyboard) mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strfind(upper(subjID),'DEBUG')

  fprintf('entering debug mode...\n');

  % just to get stimulus figures
  Screen('CloseAll');
  figure; hold on;

  save_dir=fullfile(resultDir,'images');
  if ~exist(save_dir,'dir'), mkdir(save_dir); end

  % generating and saving the depth stimuli
  imgL=cell(numel(sparam.jitters),1);
  imgR=cell(numel(sparam.jitters),1);
  for ii=1:1:length(imgL)
    imgL{ii}=cell(numel(sparam.disparity),sparam.trialsPerBlock);
    imgR{ii}=cell(numel(sparam.disparity),sparam.trialsPerBlock);
  end

  for tt=1:1:numel(sparam.jitters) % for jitter, -1,0,+1
    for ii=1:1:numel(sparam.disparity)
      for jj=1:1:sparam.trialsPerBlock
        if stim_mode==1
          [imgL{tt}{ii,jj},imgR{tt}{ii,jj}]=nf_RDSfastest(posL{tt}{ii,jj},posR{tt}{ii,jj},...
                                                          wdot,bdot,dotalpha,dotDens,sparam.colors(3),...
                                                          sparam.skipzero_flg,sparam.avoid_bias_flg,sparam.use_mex_flg);
          img=uint8([imresize(imgR{tt}{ii,jj},size(imgR{tt}{ii,jj}).*[1,1/sparam.oversampling_ratio]),...
                     sparam.colors(3)*ones(size(imgL{tt}{ii,jj},1),50),...
                     imresize(imgL{tt}{ii,jj},size(imgL{tt}{ii,jj}).*[1,1/sparam.oversampling_ratio])]);
        else
          img=CreateDepthPatchField(sparam.fieldSize,[10,10],sparam.ipd,sparam.vdist,sparam.pix_per_cm,...
                                    sparam.disparity(design(1)),[-12,12],1,sparam.oversampling_ratio,1,sparam.bgcolor);
          sz=size(img{1});
          img=uint8([imresize(img{2},sz(1:2).*[1,1/sparam.oversampling_ratio]),...
                     repmat(sparam.colors(3),[size(img{1},1),50,3]),...
                     imresize(img{1},sz(1:2).*[1,1/sparam.oversampling_ratio])]);
        end
        imshow(img); axis equal; axis off;
        fname=sprintf('wedge_depth_stimulus_jitter_%.2f_disparity_%.2f_trial_%03d.png',sparam.jitters(tt),sparam.disparity(ii),jj);
        imwrite(img,fullfile(save_dir,fname),'png');
        pause(0.1);
        close all;
      end
    end
  end
  save(fullfile(save_dir,'wedge_depth_stimulus.mat'),'imgL','imgR','sparam','dparam');
  fprintf('please check the display/stimulus parameters and the stimulus images.\n');
  keyboard;

end % if strfind(upper(subjID),'DEBUG')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Generating the first depth stimulus
%%%% this is required to present the stimulus on time at the beginning of the presentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stim=cell(2,1);
if stim_mode==1
  % wedge-shaped depth stimulus
  order=shuffle(1:1:sparam.trialsPerBlock);
  jitter_order=randi(numel(sparam.jitters),[sparam.trialsPerBlock,1]);
  [imgL,imgR]=nf_RDSfastest(posL{jitter_order(1)}{design(1),order(1)},...
                            posR{jitter_order(1)}{design(1),order(1)},...
                            wdot,bdot,dotalpha,dotDens,sparam.colors(3),...
                            sparam.skipzero_flg,sparam.avoid_bias_flg,sparam.use_mex_flg);
  stim{1}=Screen('MakeTexture',winPtr,imgL);
  stim{2}=Screen('MakeTexture',winPtr,imgR);
else
  % presenting random depth patches, just for fun,
  % ignoring the posL, posR, order, and jitter_order parameters
  % (but need to set some of those for code consistency)
  order=shuffle(1:1:sparam.trialsPerBlock);
  jitter_order=randi(numel(sparam.jitters),[sparam.trialsPerBlock,1]);
  img=CreateDepthPatchField(sparam.fieldSize,[10,10],sparam.ipd,sparam.vdist,sparam.pix_per_cm,...
                            sparam.disparity(design(1)),[-8,8],1,sparam.oversampling_ratio,1,sparam.bgcolor);
  stim{1}=Screen('MakeTexture',winPtr,img{1});
  stim{2}=Screen('MakeTexture',winPtr,img{2});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Creating the background image with vergence-guide grids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate the central aperture size of the background image
edgeY=mod(dparam.ScrHeight,sparam.patch_num(1)); % delete exceeded region
p_height=round((dparam.ScrHeight-edgeY)/sparam.patch_num(1)); % height in pix of patch_height + interval-Y

edgeX=mod(dparam.ScrWidth,sparam.patch_num(2)); % delete exceeded region
p_width=round((dparam.ScrWidth-edgeX)/sparam.patch_num(2)); % width in pix of patch_width + interval-X

aperture_size(1)=2*( p_height*ceil(size(posL{1}{1,1},1)/2/p_height) );
aperture_size(2)=2*( p_width*ceil(size(posL{1}{1,1},2)/sparam.oversampling_ratio/2/p_width) );
%aperture_size=[500,500];

bgimg=CreateBackgroundImage([dparam.ScrHeight,dparam.ScrWidth],...
          aperture_size,sparam.patch_size,sparam.bgcolor,sparam.patch_color1,sparam.patch_color2,sparam.fixcolor,sparam.patch_num,0,0,0);
background=Screen('MakeTexture',winPtr,bgimg{1});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Creating the central fixation, vernier task frame images (left/right)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create fixation cross images
[fix_L,fix_R]=CreateFixationImg(sparam.fixsize,sparam.fixcolor,sparam.bgcolor,sparam.fixlinesize(2),sparam.fixlinesize(1),0,0);
fcross{1}=Screen('MakeTexture',winPtr,fix_L);
fcross{2}=Screen('MakeTexture',winPtr,fix_R);

[wait_fix_L,wait_fix_R]=CreateFixationImg(sparam.fixsize,[32,32,32],sparam.bgcolor,sparam.fixlinesize(2),sparam.fixlinesize(1),0,0);
wait_fcross{1}=Screen('MakeTexture',winPtr,wait_fix_L);
wait_fcross{2}=Screen('MakeTexture',winPtr,wait_fix_R);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Creating task images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% creating vernier task images
vernier_texture=repmat(sparam.fixcolor(1),[sparam.verniersize,3]);

% here the *2 length of arrays are required since each presentation period is devided into 2.
vernierOffset=randsample(sparam.vernierpos,sparam.numConds*(sparam.trialsPerBlock+blank_trialsPerBlock)*sparam.numRepeats*2,true);

% to prevent image flipping problem in Haploscope settings
if strcmp(dparam.ExpMode,'Haploscope'),vernierOffset=-1*vernierOffset; end

tvernierOn=floor(rand(sparam.numConds*(sparam.trialsPerBlock+blank_trialsPerBlock)*sparam.numRepeats,1)*2);
vernierOn=[];
for ii=1:1:length(tvernierOn), vernierOn=[vernierOn,tvernierOn(ii),~tvernierOn(ii)]; end %#ok

% flags to decide in which period of each of the stimulus presentations (the first half or the second half)
% or to which (left or right) eyes, the vernier bar is presented.
vernierperiod=randi(2,[1,sparam.numConds*(sparam.trialsPerBlock+blank_trialsPerBlock)*sparam.numRepeats*2])-1;
verniereye=randi(2,[1,sparam.numConds*(sparam.trialsPerBlock+blank_trialsPerBlock)*sparam.numRepeats*2]);

% the vernier task should be presented on one of eyes
vernier=Screen('MakeTexture',winPtr,vernier_texture);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Prepare blue lines for stereo image flip sync with VPixx PROPixx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% There seems to be a blueline generation bug on some OpenGL systems.
% SetStereoBlueLineSyncParameters(winPtr, winRect(4)) corrects the
% bug on some systems, but breaks on other systems.
% We'll just disable automatic blueline, and manually draw our own bluelines!

if strcmpi(dparam.ExpMode,'propixxstereo')
  SetStereoBlueLineSyncParameters(winPtr, winRect(4)+10);
  blueRectOn(1,:)=[0, winRect(4)-1, winRect(3)/4, winRect(4)];
  blueRectOn(2,:)=[0, winRect(4)-1, winRect(3)*3/4, winRect(4)];
  blueRectOff(1,:)=[winRect(3)/4, winRect(4)-1, winRect(3), winRect(4)];
  blueRectOff(2,:)=[winRect(3)*3/4, winRect(4)-1, winRect(3), winRect(4)];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Image size adjusting to match the current display resolutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dparam.fullscr
  ratio_wid=( winRect(3)-winRect(1) )/dparam.ScrWidth;
  ratio_hei=( winRect(4)-winRect(2) )/dparam.ScrHeight;
  stimSize=[size(posL{1}{1,1},2)*ratio_wid,size(posL{1}{1,1},1)*ratio_hei].*[1/sparam.oversampling_ratio,1];
  bgSize=[size(bgimg{1},2)*ratio_wid,size(bgimg{1},1)*ratio_hei];
  fixSize=[2*sparam.fixsize*ratio_wid,2*sparam.fixsize*ratio_hei];
else
  stimSize=[size(posL{1}{1,1},2),size(posL{1}{1,1},1)].*[1/sparam.oversampling_ratio,1];
  bgSize=[dparam.ScrWidth,dparam.ScrHeight];
  fixSize=[2*sparam.fixsize,2*sparam.fixsize];
end

% for some display modes in which one screen is splitted into two binocular displays
if strcmpi(dparam.ExpMode,'cross') || strcmpi(dparam.ExpMode,'parallel') || ...
   strcmpi(dparam.ExpMode,'topbottom') || strcmpi(dparam.ExpMode,'bottomtop')
  stimSize=stimSize./2;
  bgSize=bgSize./2;
  fixSize=fixSize./2;
end

stimRect=[0,0,stimSize]; % used to display target stimuli
bgRect=[0,0,bgSize]; % used to display background images;
fixRect=[0,0,fixSize]; % used to display the central fixation point

% used to display the vernier task, as this value is so small, no adjustment depending on the resolution
vernierRect=[0,0,sparam.verniersize(2),sparam.verniersize(1)];

% set display shift along y-axis
yshift=[0,dparam.yshift,0,dparam.yshift];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Saving the current parameters temporally
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% saving the current parameters
% (this is required to analyze part of the data obtained even when the experiment is interrupted unexpectedly)
fprintf('saving the stimulus generation and presentation parameters...');
savefname=fullfile(resultDir,[num2str(subjID),sprintf('_%s_run_',mfilename()),num2str(acq,'%02d'),'.mat']);

% backup the old file(s)
if ~overwrite_flg
  rdir=relativepath(resultDir); rdir=rdir(1:end-1);
  BackUpObsoleteFiles(rdir,[num2str(subjID),sprintf('_%s_run_',mfilename()),num2str(acq,'%02d'),'.mat'],'_old');
  clear rdir;
end

% save the current parameters
eval(sprintf('save %s subjID acq design sparam dparam gamma_table;',savefname));

fprintf('done.\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Displaying 'Ready to Start'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% displaying texts on the center of the screen
%DisplayMessage2('Ready to Start',sparam.bgcolor,winPtr,nScr,'Arial',36);
ttime=GetSecs(); while (GetSecs()-ttime < 0.5), end  % run up the clock.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Flip the display(s) to the background image(s)
%%%% and inform the ready of stimulus presentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% change the screen and wait for the trigger or pressing the start button
for nn=1:1:nScr
  Screen('SelectStereoDrawBuffer',winPtr,nn-1);
  Screen('DrawTexture',winPtr,background,[],CenterRect(bgRect,winRect)+yshift);
  Screen('DrawTexture',winPtr,wait_fcross{nn},[],CenterRect(fixRect,winRect)+yshift);

  % blue line for stereo sync
  if strcmpi(dparam.ExpMode,'propixxstereo')
    Screen('FillRect',winPtr,[0,0,255],blueRectOn(nn,:));
    Screen('FillRect',winPtr,[0,0,0],blueRectOff(nn,:));
  end
end
Screen('DrawingFinished',winPtr);
Screen('Flip', winPtr,[],[],[],1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Wait for the first trigger pulse from fMRI scanner or start with button pressing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add time stamp (this also works to load add_event method in memory in advance of the actual displays)
fprintf('\nWaiting for the start...\n');
event=event.add_event('Experiment Start',strcat([datestr(now,'yymmdd'),' ',datestr(now,'HH:mm:ss')]),NaN);

% waiting for stimulus presentation
resps.wait_stimulus_presentation(dparam.start_method,dparam.custom_trigger);
%PlaySound(1);
fprintf('\nExperiment running...\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Event logs and timer (!start here!)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[event,the_experiment_start]=event.set_reference_time(GetSecs());
targetTime=the_experiment_start;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Initial Fixation Period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% wait for the initial fixation period
% it may be better to put the behavior task during this fixation period to equalize the
% attentional load with the main experiment period. But also see the comment on the task
% during the blank period below. Anyway to add tasks in the fixation period is very easy.
if sparam.initial_fixation_time(1)~=0
  event=event.add_event('Initial Fixation',[]);
  fprintf('\nfixation\n');

  for nn=1:1:nScr
    Screen('SelectStereoDrawBuffer',winPtr,nn-1);
    Screen('DrawTexture',winPtr,background,[],CenterRect(bgRect,winRect)+yshift);
    Screen('DrawTexture',winPtr,fcross{nn},[],CenterRect(fixRect,winRect)+yshift);

    % blue line for stereo sync
    if strcmpi(dparam.ExpMode,'propixxstereo')
      Screen('FillRect',winPtr,[0,0,255],blueRectOn(nn,:));
      Screen('FillRect',winPtr,[0,0,0],blueRectOff(nn,:));
    end
  end
  Screen('DrawingFinished',winPtr);
  Screen('Flip', winPtr,[],[],[],1);

  % wait for the initial fixation
  targetTime=targetTime+sparam.initial_fixation_time(1);
  while GetSecs()<targetTime, [resps,event]=resps.check_responses(event); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% The Trial Loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for currenttrial=1:1:length(design)

  cond=design(currenttrial);
  event=event.add_event('Start Block',['Cond_' int2str(cond)]);
  fprintf('Block %04d (cond %02d, %.2f arcmins): ',currenttrial,cond,sparam.disparity(cond));

  % Stimulation period
  for ii=1:sparam.trialsPerBlock

    %fprintf('% 2.1f ',sparam.jitters(jitter_order(ii)));
    fprintf('%02d ',jitter_order(ii));
    if ii==sparam.trialsPerBlock, fprintf('\n'); end

    % stimulus ON period
    % here the *2 is required since each presentation period is devided into 2.
    trialIndex=(currenttrial-1)*2*(sparam.trialsPerBlock+blank_trialsPerBlock)+2*ii-1;
    event=event.add_event('Stim On',[]);

    % devide one period into 2 for the vernier bar discrimination task
    % specifically, one period is divided into the first and the second half, and
    % the vernier task bar is presented one of the half period if vernierOn is 1.
    for jj=1:1:2
      %tStart=GetSecs();
      for nn=1:1:nScr
        Screen('SelectStereoDrawBuffer',winPtr,nn-1);
        Screen('DrawTexture',winPtr,background,[],CenterRect(bgRect,winRect)+yshift);
        Screen('DrawTexture',winPtr,stim{nn},[],CenterRect(stimRect,winRect)+yshift);
        Screen('DrawTexture',winPtr,fcross{nn},[],CenterRect(fixRect,winRect)+yshift);

        % vernier task
        if vernierOn(trialIndex+jj-1) && vernierperiod(trialIndex+jj-1) && verniereye(trialIndex+jj-1)==nn
          event=event.add_event('Vernier',vernierOffset(trialIndex+jj-1));
          Screen('DrawTexture',winPtr,vernier,[],...
                 CenterRect(vernierRect,winRect)+[vernierOffset(trialIndex),0,vernierOffset(trialIndex),0]+yshift);
        end

        % blue line for stereo sync
        if strcmpi(dparam.ExpMode,'propixxstereo')
          Screen('FillRect',winPtr,[0,0,255],blueRectOn(nn,:));
          Screen('FillRect',winPtr,[0,0,0],blueRectOff(nn,:));
        end
      end
      Screen('DrawingFinished',winPtr);
      Screen('Flip',winPtr,[],[],[],1);

      % wait for stim_on_duration
      % here the /2 is required since each presentation period is devided into 2.
      targetTime=targetTime+sparam.stim_on_duration/2;
      while GetSecs()<targetTime, [resps,event]=resps.check_responses(event); end
    end % for jj=1:1:2

    % clean up the current texture & release memory
    for mm=1:1:2, Screen('Close',stim{mm}); end

    % back to the background only display
    for nn=1:1:nScr
      Screen('SelectStereoDrawBuffer',winPtr,nn-1);
      Screen('DrawTexture',winPtr,background,[],CenterRect(bgRect,winRect)+yshift);
      Screen('DrawTexture',winPtr,fcross{nn},[],CenterRect(fixRect,winRect)+yshift);

      % blue line for stereo sync
      if strcmpi(dparam.ExpMode,'propixxstereo')
        Screen('FillRect',winPtr,[0,0,255],blueRectOn(nn,:));
        Screen('FillRect',winPtr,[0,0,0],blueRectOff(nn,:));
      end
    end
    Screen('DrawingFinished',winPtr);
    Screen('Flip',winPtr,[],[],[],1);

    % preparing the next stimuli

    % check & update the next trial index.
    % if it comes to the final condition, only background is presented for stim_off_duration
    if ii<sparam.trialsPerBlock % stay in the current condition
      nexttrial=currenttrial;
      nextcond=design(nexttrial);
    else % go to the next condition
      nexttrial=currenttrial+1;
      % wait for the final fixation period
      if nexttrial>length(design)
        targetTime=targetTime +sparam.stim_off_duration;
        while GetSecs()<targetTime, [resps,event]=resps.check_responses(event); end
        break;
      else
        nextcond=design(nexttrial);
      end
    end

    [resps,event]=resps.check_responses(event); % just in case, get the response also here

    % generate next RDS
    if stim_mode==1
      % wedge-shaped depth stimulus
      if ii<sparam.trialsPerBlock % stay in the current condition
        [imgL,imgR]=nf_RDSfastest(posL{jitter_order(ii+1)}{nextcond,order(ii+1)},...
                                  posR{jitter_order(ii+1)}{nextcond,order(ii+1)},...
                                  wdot,bdot,dotalpha,dotDens,sparam.colors(3),...
                                  sparam.skipzero_flg,sparam.avoid_bias_flg,sparam.use_mex_flg);
      else
        % generate the next order array to randomize the sequence of posL/R{disparity,trialsPerBlock}
        order=shuffle(1:1:sparam.trialsPerBlock);
        jitter_order=randi(numel(sparam.jitters),[sparam.trialsPerBlock,1]);
        [imgL,imgR]=nf_RDSfastest(posL{jitter_order(1)}{nextcond,order(1)},...
                                  posR{jitter_order(1)}{nextcond,order(1)},...
                                  wdot,bdot,dotalpha,dotDens,sparam.colors(3),...
                                  sparam.skipzero_flg,sparam.avoid_bias_flg,sparam.use_mex_flg);
      end
      stim{1}=Screen('MakeTexture',winPtr,imgL);
      stim{2}=Screen('MakeTexture',winPtr,imgR);
    else
      % presenting random depth patches, just for fun,
      % ignoring the posL, posR, order, and jitter_order parameters
      % (but need to set some of those for code consistency)
      order=shuffle(1:1:sparam.trialsPerBlock);
      jitter_order=randi(numel(sparam.jitters),[sparam.trialsPerBlock,1]);
      img=CreateDepthPatchField(sparam.fieldSize,[10,10],sparam.ipd,sparam.vdist,sparam.pix_per_cm,...
                                sparam.disparity(nextcond),[-8,8],1,sparam.oversampling_ratio,1,sparam.bgcolor);
      stim{1}=Screen('MakeTexture',winPtr,img{1});
      stim{2}=Screen('MakeTexture',winPtr,img{2});
    end

    % wait for stim_off_duration
    targetTime=targetTime+sparam.stim_off_duration;
    while GetSecs()<targetTime, [resps,event]=resps.check_responses(event); end

  end % for ii=1:dparam.trialsPerBlock

  % Blank period
  %
  % note:
  % here, generally speaking, the same vernier task performed during the stimulation period
  % may be required even during the blank period. However, the vernier task is used to check
  % whether the participant vergences are fine when the target stimuli are presented and
  % whether the participant vergences are consistent across the stimulus conditions. To those
  % purposes, what we want to know is accuracies of the task performed during stimulus
  % presentation periods only. Rather, if we add the vernier task accuracies during the
  % blank period, it may cause some problems since the perfromance may be obscured as the
  % overall accuracies would be boosted during the blank period in which participants can
  % concentrate on the tasks more since no depth stimulus is provided in the period.
  % Therefore, I intentionally omit the task during the blank period. This decision may
  % cause attention load differences between the stimulation and blank period. Please
  % be sure this is fine for your experiment.
  %
  % if you want to add the venier task during the blank period, please uncomment the
  % vernier task codes below.

  if sparam.blank_duration~=0
    event=event.add_event('Blank',[]);
    fprintf('Blank %04d (cond %02d, no stimulus): ... ',currenttrial,cond);
    for ii=1:1:blank_trialsPerBlock

      % stimulus ON period
      % here the *2 is required since each presentation period is devided into 2.
      %trialIndex=(currenttrial-1)*2*(sparam.trialsPerBlock+blank_trialsPerBlock)+sparam.trialsPerBlock+2*ii-1;

      for jj=1:1:2 % devide one period into 2 for the vernier bar discrimination task
        %tStart=GetSecs();
        for nn=1:1:nScr
          Screen('SelectStereoDrawBuffer',winPtr,nn-1);
          Screen('DrawTexture',winPtr,background,[],CenterRect(bgRect,winRect)+yshift);
          Screen('DrawTexture',winPtr,fcross{nn},[],CenterRect(fixRect,winRect)+yshift);

          % *** intentionally disabled now ***
          % % vernier task
          % if vernierOn(trialIndex+jj-1) && vernierperiod(trialIndex+jj-1) && verniereye(trialIndex+jj-1)==nn
          %   event=event.add_event('Vernier',vernierOffset(trialIndex+jj-1));
          %   Screen('DrawTexture',winPtr,vernier,[],...
          %          CenterRect(vernierRect,winRect)+[vernierOffset(trialIndex),0,vernierOffset(trialIndex),0]+yshift);
          % end

          % blue line for stereo sync
          if strcmpi(dparam.ExpMode,'propixxstereo')
            Screen('FillRect',winPtr,[0,0,255],blueRectOn(nn,:));
            Screen('FillRect',winPtr,[0,0,0],blueRectOff(nn,:));
          end
        end
        Screen('DrawingFinished',winPtr);
        Screen('Flip',winPtr,[],[],[],1);

        % wait for stim_on_duration
        % here the /2 is required since each presentation period is devided into 2.
        targetTime=targetTime+sparam.stim_on_duration/2;
        while GetSecs()<targetTime, [resps,event]=resps.check_responses(event); end
      end % for jj=1:1:2

      % wait for stim_off_duration
      targetTime=targetTime+sparam.stim_off_duration;
      while GetSecs()<targetTime, [resps,event]=resps.check_responses(event); end

    end % for ii=1:1:blank_trialsPerBlock
    fprintf('\n');
  end % if sparam.blank_duration~=0

end % for currenttrial=1:1:length(design)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Final Fixation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% wait for the initial fixation
% it may be better to put the behavior task during this fixation period to equalize the
% attentional load with the main experiment period. But also see the comment on the task
% during the blank period above. Anyway to add tasks in the fixation period is very easy.
if sparam.initial_fixation_time(2)~=0
  event=event.add_event('Final Fixation',[]);
  fprintf('fixation\n');

  for nn=1:1:nScr
    Screen('SelectStereoDrawBuffer',winPtr,nn-1);
    Screen('DrawTexture',winPtr,background,[],CenterRect(bgRect,winRect)+yshift);
    Screen('DrawTexture',winPtr,fcross{nn},[],CenterRect(fixRect,winRect)+yshift);

    % blue line for stereo sync
    if strcmpi(dparam.ExpMode,'propixxstereo')
      Screen('FillRect',winPtr,[0,0,255],blueRectOn(nn,:));
      Screen('FillRect',winPtr,[0,0,0],blueRectOff(nn,:));
    end
  end
  Screen('DrawingFinished',winPtr);
  Screen('Flip', winPtr,[],[],[],1);

  % wait for the initial fixation
  targetTime=targetTime+sparam.initial_fixation_time(2);
  while (GetSecs()<targetTime), [resps,event]=resps.check_responses(event); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Experiment & scanner end here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

experimentDuration=GetSecs()-the_experiment_start;
event=event.add_event('End',[]);

fprintf('\nExperiment Completed: %.2f/%.2f secs\n\n',experimentDuration,...
  sum(sparam.initial_fixation_time)+sparam.numConds*sparam.numRepeats*sparam.block_duration);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Write data into a file for further analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% saving the results
fprintf('saving data...');

% save data
savefname=fullfile(resultDir,[num2str(subjID),sprintf('_%s_run_',mfilename()),num2str(acq,'%02d'),'.mat']);
eval(sprintf('save -append %s subjID acq sparam dparam design event gamma_table;',savefname));
fprintf('done.\n');

% calculate & display task performance
% The codes below simply calculates the overall task accuracies, ignoring the condition differnces.
% if you want to calculate the accuracies by your own criteria etc, please change the codes below.
correct_events=cell(numel(sparam.vernierpos),1);
for ii=1:1:numel(sparam.vernierpos)
  if sparam.vernierpos(ii)<=0
    correct_events{ii}={sparam.vernierpos(ii),'key1'};
  else
    correct_events{ii}={sparam.vernierpos(ii),'key2'};
  end
end
[task.numTasks,task.numHits,task.numErrors,task.numResponses,task.RT,event]=event.calc_accuracy(correct_events);

eval(sprintf('save -append %s event task;',savefname));
fprintf('done.\n');

% tell the experimenter that the measurements are completed
try
  for ii=1:1:3, Snd('Play',sin(2*pi*0.2*(0:900)),8000); end
catch
  % do nothing
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Cleaning up the PTB screen, removing path to the subfunctions, and finalizing the script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Screen('CloseAll');

% closing datapixx
if strcmpi(dparam.exp_mode,'propixxmono') || strcmpi(dparam.exp_mode,'propixxstereo')
  if Datapixx('IsViewpixx3D')
    Datapixx('DisableVideoLcd3D60Hz');
    Datapixx('RegWr');
  end
  Datapixx('Close');
end

ShowCursor();
Priority(0);
GammaResetPTB(1.0);
rmpath(genpath(fullfile(rootDir,'..','Common')));
rmpath(fullfile(rootDir,'..','Generation'));
clear all; clear mex; clear global;
diary off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Catch the errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

catch %#ok
  % this "catch" section executes in case of an error in the "try" section
  % above.  Importantly, it closes the onscreen window if its open.
  Screen('CloseAll');

  if exist('dparam','var')
    if isstructmember(dparam,'exp_mode')
      if strcmpi(dparam.exp_mode,'propixxmono') || strcmpi(dparam.exp_mode,'propixxstereo')
        if Datapixx('IsViewpixx3D')
          Datapixx('DisableVideoLcd3D60Hz');
          Datapixx('RegWr');
        end
        Datapixx('Close');
      end
    end
  end

  ShowCursor();
  Priority(0);
  GammaResetPTB(1.0);
  tmp=lasterror; %#ok
  if exist('event','var'), event=event.get_event(); end %#ok % just for debugging
  diary off;
  fprintf(['\nError detected and the program was terminated.\n',...
           'To check error(s), please type ''tmp''.\n',...
           'Please save the current variables now if you need.\n',...
           'Then, quit by ''dbquit''\n']);
  keyboard;
  rmpath(genpath(fullfile(rootDir,'..','Common')));
  rmpath(fullfile(rootDir,'..','Generation'));
  clear all; clear mex; clear global;
  return
end % try..catch


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% That's it - we're done
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return
% end % function StereofMRIsample
