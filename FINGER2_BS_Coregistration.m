%% Coregistration of between-session group in Arana et al. (2025).

% This script is employed for semi-automatic coregistration of data from The Open MEG Archive (Niso et al., 2016). 
% The output is the mri_coreg.
% The input is the T1 mri of each subject and their headshapes.
% You should select the first session or the second.

clear all
close all
clc

restoredefaultpath
addpath ('Z:\Toolbox\fieldtrip-20230118');
ft_defaults
addpath(genpath('Z:\Fingerprinting\scripts\Final'));


dpath = 'Z:\OMEGA\OMEGA_data\';
rawpath = 'Z:\OMEGA\OMEGA_raw\';

% Subjects of between-session group that have 2 sessions
subs = {'sub-0001'	'sub-0002'	'sub-0006'	'sub-0008'	'sub-0011'	'sub-0016'	'sub-0019'	'sub-0020'	'sub-0022'	'sub-0023'	'sub-0025'	'sub-0030'	'sub-0032'	'sub-0035'	'sub-0039'	'sub-0040'	'sub-0041'	'sub-0042'	'sub-0044'	'sub-0046'	'sub-0048'	'sub-0049'	'sub-0050'	'sub-0051'	'sub-0106'	'sub-0150'	'sub-0200'};

% Select session
select_session = 1;   % 1 for fisrt session, 2 for second session

if select_session == 1

    % First sess
    sess = {'ses-0001' 'ses-0004' 'ses-0001'  'ses-0001'	'ses-0002'	'ses-0002'	'ses-0003'	'ses-0004'	'ses-0002'	'ses-0001'	'ses-0003'	'ses-0001'	'ses-0001'	'ses-0003'	'ses-0001'	'ses-0005'	'ses-0003'	'ses-0001'  'ses-0002'	'ses-0002'	'ses-0004'	'ses-0002'	'ses-0004'	'ses-0001'  'ses-0001'	'ses-0001'	'ses-0001'};

elseif select_session == 2

    % Second sess
    sess = {'ses-0003' 'ses-0002' 'ses-0003'	'ses-0002'	'ses-0001' 	'ses-0001'	'ses-0002'	'ses-0003'	'ses-0001'	'ses-0002'	'ses-0002'	'ses-0002'	'ses-0002'	'ses-0002'	'ses-0002'	'ses-0004'	'ses-0001'	'ses-0002'	'ses-0003'	'ses-0001'	'ses-0002'	'ses-0001'	'ses-0002'	'ses-0002'	'ses-0002'	'ses-0002'	'ses-0002'};

else
    disp('You should select session 1 or 2')
end


% 2.1. Coregistration of MEG-MRI spaces (output: mri_coreg)

for sub = 1:length(subs)
%% 2.1. Coregistration of MEG-MRI spaces
% mri.transform is the transformation matrix to go from mri space to sensor space

cd([dpath subs{sub}])
mri = ft_read_mri('defaced_t1.nii');

cd([rawpath  subs{sub} '\' sess{sub} '\meg'])

dataresting = findfile('resting');
cd(dataresting)                                     

try
    hsfile    = findfile('.pos');
    headshape = ft_read_headshape(hsfile);
catch
    hsfile    = findfile('corrected.pos');      % corrected hs files (no letters in 2nd column)  
    headshape = ft_read_headshape(hsfile);
end
[mri,scp] = omega_coreg([], mri, headshape);    % mark fiducials: lpa (l), rpa (r) and nasion (n), then quit (q) 

cd([dpath  subs{sub} '\' sess{sub} ])
savefig(gcf, 'corr.fig');
save mri_coreg mri scp
end

