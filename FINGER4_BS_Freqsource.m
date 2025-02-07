%% Freqsource of between-session group in Arana et al. (2025).

% This script is employed for obtention of power spectra from continuous data of The Open MEG Archive. 
% The outputs are the power spectra of each session freq_allvox_10mm_steps100ms.
% The inputs are dataclean, source_forward_10mm and source_inverse_10mm.

clear all
close all
clc

restoredefaultpath
addpath ('Z:\Toolbox\fieldtrip-20230118');
ft_defaults
addpath(genpath('Z:\Fingerprinting\scripts\Final'));

dpath = 'Z:\OMEGA\OMEGA_data\';
rawpath = 'Z:\OMEGA\OMEGA_raw\';



subs = {'sub-0001'	'sub-0002'	'sub-0006'	'sub-0008'	'sub-0011'	'sub-0016'	'sub-0019'	'sub-0020'	'sub-0022'	'sub-0023'	'sub-0025'	'sub-0030'	'sub-0032'	'sub-0035'	'sub-0039'	'sub-0040'	'sub-0041'	'sub-0042'	'sub-0044'	'sub-0046'	'sub-0048'	'sub-0049'	'sub-0050'	'sub-0051'	'sub-0106'	'sub-0150'	'sub-0200'};

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


% 4.1. Reconstruction of source-level time-series
% 4.2. Frequency analysis parameters
% 4.3. Frequency analysis computation (output: freq_allvox_10mm_steps100ms)

for sub=1:length(subs)
%% 4.1. Reconstruction of source-level time-series

% Read necessary files
cd([dpath  subs{sub} '\' sess{sub} ])
load dataclean
load source_forward_10mm
load source_inverse_10mm

% Reconstruct source-space data
time         = dataclean.time{1};
voxel_inside = find(source.inside==1);
Nvox         = length(voxel_inside);
datasource   = zeros(Nvox,length(time));
for i = 1:Nvox
    disp([num2str(i) ' / ' num2str(length(voxel_inside))])
    datasource(i,:) = source.avg.filter{voxel_inside(i)} * dataclean.trial{1};
end

datasource = datasource./repmat(std(datasource,0,2),[1,length(time)]);        % baseline correction to account for the centre of the head bias

% Discard artifactual time segments 
cd([dpath  subs{sub} '\' sess{sub} ])
if exist('badsegments.mat') == 2
    load badsegments
    for b = 1:length(badsegments)
        t1 = findbin(time,badsegments{b}(1));
        t2 = findbin(time,badsegments{b}(2));
        time(t1:t2) = [];
        datasource(:,t1:t2) = [];
    end
end

% Use only 5-minute recordings for all participants
if time(end) > 290      % >5 min
    t1 = findbin(time,290);
    t2 = findbin(time,time(end));
    time(t1:t2) = [];
    datasource(:,t1:t2) = [];
end

% Organize source-reconstructed data in a new fieldtrip structure
dataclean.trial{1}   = single(datasource);
dataclean.time{1}    = time;
dataclean.label      = {};
dataclean.sampleinfo = [1 size(dataclean.trial{1},2)];
for i = 1:Nvox
    dataclean.label{i} = ['V' num2str(i)];
end

clear datasource

% If there were not badsegments, data is organized in 1 single trial
% If there were badsegments (discontinuities), data is organized in several trials
if exist('badsegments.mat')==2
    dtime      = diff(time);
    [pks,locs] = findpeaks(dtime);
    pkstime    = time(locs);
    if length(pkstime) > 0
        trl = [];
        trl(1,1) = 1;
        for t = 1:length(pks)
            tt = findbin(time, pkstime(t));
            trl(t,2) = tt - 1;
            trl(t+1,1) = tt + 1;
        end
        trl(t+1,2) = length(time);
        trl(:,3)   = 0;
        
        cfg     = [];
        cfg.trl = trl;
        dataclean = ft_redefinetrial(cfg,dataclean);     
    end
end


%% 4.2. Frequency analysis parameters: foi, toi, t_ftimwin

f         = 0.55:0.05:3.55;        % originally up to 4.6 (now up to 34.8 Hz)
foi       = exp(f);                % logarithmically spaced frequencies
t_ftimwin = 5./foi;                % time-window adapted to each frequency (5 cycles)
t_fstep   = 0.1;                   % sliding time-window (100 ms steps)

bd   = t_ftimwin(1).*2;            % remove borders = 2*time-window
bdpt = bd.*dataclean.fsample;

Ntr  = 0;
toi2 = {};
ct   = 1;
valid_tr = [];
for tr = 1:length(dataclean.trial)
    toi = bd:t_fstep:dataclean.time{tr}(end)-bd;
    if length(toi) > 0
        toi2{ct} = toi;
        Ntr = Ntr+length(toi);
        valid_tr(ct) = tr;
        ct = ct+1;
    end
end


%% 4.3. Frequency analysis computation

Nvox  = length(dataclean.label);
Nf    = length(foi);
powsp = single(NaN(Nvox,Nf,Ntr));

tr2 = 0;
for tr = 1:length(valid_tr)
    tr1 = tr2+1;
    tr2 = tr1+length(toi2{tr})-1;
    
    cfg            = [];
    cfg.trials     = valid_tr(tr);
    cfg.method     = 'mtmconvol';
    cfg.taper      = 'hanning';
    cfg.output     = 'pow';
    cfg.foi        = foi;
    cfg.toi        = toi2{tr};
    cfg.t_ftimwin  = t_ftimwin;
    cfg.pad        = 'nextpow2';
    cfg.keeptrials = 'yes';
    freq           = ft_freqanalysis(cfg, dataclean);
    
    powsp(:,:,tr1:tr2)= single(squeeze(freq.powspctrm));
end

trok=[];
ct=1;
for tr = 1:size(powsp,3)
    if sum(sum(isnan(powsp(:,:,tr))))==0
        trok(ct) = tr;
        ct=ct+1;
    end
end

powsp = powsp(:,:,trok);
Ntr_ok=length(trok) % all trials of trials/spectra of the subj

cd([dpath  subs{sub} '\' sess{sub} ])
save freq_allvox_10mm_steps100ms powsp foi Ntr_ok

end

