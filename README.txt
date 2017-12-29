**************************************************
README.txt on StereofMRIsample

Created    : "2017-12-29 13:08:04 ban"
Last Update: "2017-12-29 14:58:56 ban"
**************************************************


A block-design fMRI experiment with wedge-shaped stimulu consisted of near/far random-dot-stereograms (RDSs).
function StereofMRIsample(subjID,acq,:displayfile,:stimulusfile,:gamma_table,:overwrite_flg)
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

For details, please see the documents in StereofMRIsample.m
Also please see the parameter files in ~/StereofMRIsample/Presentation/subj/TEST.
