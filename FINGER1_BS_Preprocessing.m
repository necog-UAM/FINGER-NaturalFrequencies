%% Preprocesssing of between-session group in Arana et al. (2025).

% This script is employed for manually preprocessing data from The Open MEG Archive (Niso et al., 2016).
% The main outputs are the clean continuous recording after ICA removal (dataclean) and badsegments saving (badsegments) for avoiding them later.
% The input to use this script is the raw data from The Open MEG Archive.
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

    % First session
    sess = {'ses-0001' 'ses-0004' 'ses-0001'  'ses-0001'	'ses-0002'	'ses-0002'	'ses-0003'	'ses-0004'	'ses-0002'	'ses-0001'	'ses-0003'	'ses-0001'	'ses-0001'	'ses-0003'	'ses-0001'	'ses-0005'	'ses-0003'	'ses-0001'  'ses-0002'	'ses-0002'	'ses-0004'	'ses-0002'	'ses-0004'	'ses-0001'  'ses-0001'	'ses-0001'	'ses-0001'};

elseif select_session == 2

    % Second session
    sess = {'ses-0003' 'ses-0002' 'ses-0003'	'ses-0002'	'ses-0001' 	'ses-0001'	'ses-0002'	'ses-0003'	'ses-0001'	'ses-0002'	'ses-0002'	'ses-0002'	'ses-0002'	'ses-0002'	'ses-0002'	'ses-0004'	'ses-0001'	'ses-0002'	'ses-0003'	'ses-0001'	'ses-0002'	'ses-0001'	'ses-0002'	'ses-0002'	'ses-0002'	'ses-0002'	'ses-0002'};

else
    disp('You should select session 1 or 2')
end


% 1.1. Reading data
% 1.2. Preprocessing (output: data)
% 1.3. Artifact correction with ICA (output: dataclean, comp, badsegments)

for sub = 1:length(subs)

%% 1.1. Reading data

cd([rawpath  subs{sub} '\' sess{sub} '\meg'])

dataresting = findfile('resting');
cd(dataresting)                                  

datafile    = findfile('.meg4');
headerfile  = findfile('.res4');

cfg            = [];
cfg.datafile   = datafile;
cfg.headerfile = headerfile;
cfg.continuous = 'yes';
data           = ft_preprocessing(cfg);

%% 1.2. Preprocessing 

cfg            = [];
cfg.refchannel = 'MEGREF';
data           = ft_denoise_pca(cfg,data);

cfg            = [];
cfg.hpfilter   = 'yes';
cfg.hpfreq     = 0.05;
cfg.hpfiltord  = 3;
cfg.dftfilter  = 'yes';
cfg.dftfreq    = [50 100 150];
cfg.dftreplace = 'neighbour_fft';
data           = ft_preprocessing(cfg,data);

cfg            = [];
cfg.resamplefs = 512;
cfg.detrend    = 'yes';
cfg.demean     = 'yes';
data           = ft_resampledata(cfg,data);

cfg         = [];
cfg.toilim  = [10 data.time{1}(end)-10];
data        = ft_redefinetrial(cfg,data);

cd([dpath  subs{sub} '\' ])
mkdir([dpath  subs{sub} '\' sess{sub} ])
cd([dpath  subs{sub} '\' sess{sub} ])
save data data

%% 1.3. Artifact correction with ICA 

cfg              = [];
cfg.numcomponent = 40;                          % keep 40 principal components               
cfg.method       = 'runica';
comp             = ft_componentanalysis(cfg ,data);         

cd([dpath  subs{sub} '\' sess{sub} ])
save comp comp

% Detect badicas and identify badsegments (many ICs contaminated during the same time segment)
cfg        = [];
cfg.layout = 'CTF275.lay';
cfg.scale  = 4e-11; 
badicas    = detect_badicas_omega (cfg, data, comp);      % automatically save badicas: mouse-click reject; space-bar keep component
% badicas = [];
% save badicas badicas                                    % manually save badicas 

% Remove badicas from signal and save dataclean
cfg           = [];
cfg.component = badicas;
dataclean     = ft_rejectcomponent (cfg, comp, data);

cfg             = [];
cfg.viewmode    = 'vertical'; 
cfg.ylim = [-6.25e-13    6.25e-13];
cfg.blocksize=12;
ft_databrowser (cfg, dataclean); % Visualize segments

% figure, plot(dataclean.time{1},dataclean.trial{1})    % example for visual inspection of data to find badsegments
% ica=15
% figure, plot(comp.time{1},comp.trial{1}(ica,:))       % example for visualize components to find badsegments

% badsegments{1}= [183 187]                             % example for saving badsegments 
% badsegments{2}= [241 246] 
% badsegments{3}= [198 199]  
% badsegments{4}= [133 138] 
% badsegments{5}= [237 240] 
% cd([dpath  subs{sub} '\' sess{sub} ])
% save badsegments badsegments 

cd([dpath  subs{sub} '\' sess{sub} ])
save dataclean dataclean
end
