function OK=run_exp(subj,exp_num,run_ids)

% function run_exp(subj,exp_num,run_ids)
% (: is optional)
%
% run StereofMRIsample -- near/far RDSs (Random Dot-Stereograms)
%
% [input]
% subj    : subject's name, e.g. 'HB'
% exp_num : experiment ID that you want to run,
%           should be 1 (single zero-order disparity wedge) or
%                     2 (planes with random disparities)
% run_ids : ids you want to run, 1,2,3,...
%           a scalar or a [1xN] matrix.
%
% [output]
% OK      : whether this script is performed
%           without any error [true/false]
%
% Created    : "2017-12-29 13:07:06 ban"
% Last Update: "2019-03-01 15:15:30 ban"

%% input variable check
if nargin<3, help(mfilenae()); return; end

if isempty(intersect(exp_num,[1,2]))
  warning('MATLAB:exp_num_error','exp_num should be 1 or 2. please chech the input variable.');
  OK=false;
  return;
end

%% set stimulus file name
%if exp_num<=2
%  stim_fname=sprintf('stimulus_exp%d',exp_num);
%else
%  stim_fname='nearfar_stimulus_mixed';
%end

%% check directory with subject name

% [NOTE]
% if the subj directory is not found, create subj directory, copy all condition files
% from _DEFAULT_, and then run the script using the parameters in the default directory
subj_dir=fullfile(fileparts(mfilename('fullpath')),'subjects',subj);
if ~exist(subj_dir,'dir')
  disp('The subject directory was not found.');
  user_response=0;
  while ~user_response
    user_entry = input('Do you want to proceed using DEFAULT parameters? (y/n) : ', 's');
    if(user_entry == 'y')
      fprintf('Generating subj directory using DEFAULT parameters...');
      user_response=1; %#ok
      break;
    elseif (user_entry == 'n')
      disp('quiting the script...');
      if nargout, OK=false; end
      return;
    else
      disp('Please answer y or n!'); continue;
    end
  end

  %mkdir(subj_dir);
  copyfile(fullfile(fileparts(mfilename('fullpath')),'subjects','_DEFAULT_'),subj_dir);
end

% ********************************************************************************************************
% *** set gamma table. please change the line below to use the actual measuments of the display gamma. ***
% ********************************************************************************************************

% loading gamma_table
load(fullfile('..','gamma_table','ASUS_ROG_Swift_PG278Q','181003','cbs','gammatablePTB.mat'));
%load(fullfile('..','gamma_table','ASUS_VG278HE','181003','cbs','gammatablePTB.mat'));
%load(fullfile('..','gamma_table','MEG_B1','151225','cbs','gammatablePTB.mat'));
%gammatable=repmat(linspace(0.0,1.0,256),3,1)'; %#ok % a simple linear gamma

%% run the stimulus presentation code
for ii=run_ids
  command_str=sprintf('StereofMRIsample(''%s'',%d,''sample_displayfile'',''sample_stimulusfile'',gammatable,1,%d);',subj,ii,exp_num);
  eval(command_str);
end
OK=true;

return
