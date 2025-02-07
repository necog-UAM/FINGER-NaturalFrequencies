%% Preprocesssing of within-session group in Arana et al. (2025).

% This script is employed for manually preprocessing data from The Open MEG Archive (Niso et al., 2016). 
% The main outputs are the clean continuous recording after ICA removal (dataclean) and badsegments saving (badsegments) for avoiding them later.
% The input to use this script is the raw data from The Open MEG Archive.

clear all
close all
clc

restoredefaultpath
addpath ('Z:\Toolbox\fieldtrip-20230118');
ft_defaults
addpath(genpath('Z:\Fingerprinting\scripts\Final'));


dpath = 'Z:\OMEGA\OMEGA_data\';
rawpath = 'Z:\OMEGA\OMEGA_raw\';


% One session, within-session group
subs = {'sub-0001' 'sub-0002' 'sub-0003' 'sub-0004' 'sub-0005' 'sub-0006' 'sub-0007' 'sub-0008' 'sub-0009' 'sub-0011' 'sub-0012' 'sub-0014' 'sub-0015' 'sub-0016' 'sub-0018' 'sub-0019' 'sub-0020' 'sub-0021' 'sub-0022' 'sub-0023' 'sub-0024' 'sub-0025' 'sub-0026' 'sub-0027' 'sub-0028' 'sub-0029' 'sub-0030' 'sub-0031' 'sub-0032' 'sub-0033' 'sub-0034' 'sub-0035' 'sub-0037' 'sub-0039' 'sub-0040' 'sub-0041' 'sub-0042' 'sub-0044' 'sub-0045' 'sub-0046' 'sub-0047' 'sub-0048' 'sub-0049' 'sub-0050' 'sub-0051' 'sub-0052' 'sub-0055' 'sub-0056' 'sub-0057' 'sub-0058' 'sub-0059' 'sub-0060' 'sub-0061' 'sub-0062' 'sub-0063' 'sub-0064' 'sub-0065' 'sub-0067' 'sub-0068' 'sub-0069' 'sub-0070' 'sub-0071' 'sub-0072' 'sub-0073' 'sub-0074' 'sub-0075' 'sub-0076' 'sub-0077' 'sub-0078' 'sub-0079' 'sub-0080' 'sub-0084' 'sub-0085' 'sub-0087' 'sub-0088' 'sub-0089' 'sub-0090' 'sub-0091' 'sub-0092' 'sub-0094' 'sub-0095' 'sub-0096' 'sub-0097' 'sub-0098' 'sub-0099' 'sub-0101' 'sub-0102' 'sub-0103' 'sub-0104' 'sub-0105' 'sub-0106' 'sub-0134' 'sub-0145' 'sub-0146' 'sub-0148' 'sub-0149' 'sub-0150' 'sub-0151' 'sub-0152' 'sub-0154' 'sub-0155' 'sub-0156' 'sub-0157' 'sub-0158' 'sub-0159' 'sub-0160' 'sub-0161' 'sub-0165' 'sub-0166' 'sub-0167' 'sub-0168' 'sub-0169' 'sub-0170' 'sub-0171' 'sub-0175' 'sub-0176' 'sub-0177' 'sub-0179' 'sub-0181' 'sub-0184' 'sub-0185' 'sub-0195' 'sub-0197' 'sub-0200' 'sub-0207' 'sub-0208' 'sub-0210' 'sub-0212'};
sess = {'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0002' 'ses-0003' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0003' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001'};


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
