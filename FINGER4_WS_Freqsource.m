%% Freqsource of within-session group in Arana et al. (2025).

% This script is employed for obtention of power spectra from continuous data of The Open MEG Archive. 
% Specifically, this was employed for frequency analysis of within-session group that uses the both halves of one MEG session for each subject, separately. 
% The outputs are the power spectra of each half freq_allvox_10mm_steps100ms_1 and freq_allvox_10mm_steps100ms_2.
% The inputs are dataclean1, source_forward_10mm_1 and source_inverse_10mm_1 
% or dataclean2, source_forward_10mm_2 and source_inverse_10mm_2. 

clear all
close all
clc

restoredefaultpath
addpath ('Z:\Toolbox\fieldtrip-20230118');
ft_defaults
addpath(genpath('Z:\Fingerprinting\scripts\Final'));


dpath = 'Z:\OMEGA\OMEGA_data\';
outpath = 'G:\Fingerprinting\Omega_data\';


% One sess, within-session group
subs = {'sub-0001' 'sub-0002' 'sub-0003' 'sub-0004' 'sub-0005' 'sub-0006' 'sub-0007' 'sub-0008' 'sub-0009' 'sub-0011' 'sub-0012' 'sub-0014' 'sub-0015' 'sub-0016' 'sub-0018' 'sub-0019' 'sub-0020' 'sub-0021' 'sub-0022' 'sub-0023' 'sub-0024' 'sub-0025' 'sub-0026' 'sub-0027' 'sub-0028' 'sub-0029' 'sub-0030' 'sub-0031' 'sub-0032' 'sub-0033' 'sub-0034' 'sub-0035' 'sub-0037' 'sub-0039' 'sub-0040' 'sub-0041' 'sub-0042' 'sub-0044' 'sub-0045' 'sub-0046' 'sub-0047' 'sub-0048' 'sub-0049' 'sub-0050' 'sub-0051' 'sub-0052' 'sub-0055' 'sub-0056' 'sub-0057' 'sub-0058' 'sub-0059' 'sub-0060' 'sub-0061' 'sub-0062' 'sub-0063' 'sub-0064' 'sub-0065' 'sub-0067' 'sub-0068' 'sub-0069' 'sub-0070' 'sub-0071' 'sub-0072' 'sub-0073' 'sub-0074' 'sub-0075' 'sub-0076' 'sub-0077' 'sub-0078' 'sub-0079' 'sub-0080' 'sub-0084' 'sub-0085' 'sub-0087' 'sub-0088' 'sub-0089' 'sub-0090' 'sub-0091' 'sub-0092' 'sub-0094' 'sub-0095' 'sub-0096' 'sub-0097' 'sub-0098' 'sub-0099' 'sub-0101' 'sub-0102' 'sub-0103' 'sub-0104' 'sub-0105' 'sub-0106' 'sub-0134' 'sub-0145' 'sub-0146' 'sub-0148' 'sub-0149' 'sub-0150' 'sub-0151' 'sub-0152' 'sub-0154' 'sub-0155' 'sub-0156' 'sub-0157' 'sub-0158' 'sub-0159' 'sub-0160' 'sub-0161' 'sub-0165' 'sub-0166' 'sub-0167' 'sub-0168' 'sub-0169' 'sub-0170' 'sub-0171' 'sub-0175' 'sub-0176' 'sub-0177' 'sub-0179' 'sub-0181' 'sub-0184' 'sub-0185' 'sub-0195' 'sub-0197' 'sub-0200' 'sub-0207' 'sub-0208' 'sub-0210' 'sub-0212'};
sess = {'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0002' 'ses-0003' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0003' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001'};


% 4.1. Reconstruction of source-level time-series
% 4.2. Frequency analysis parameters
% 4.3. Frequency analysis computation (output: freq_allvox_10mm_steps100ms_1 and freq_allvox_10mm_steps100ms_2 )

for sub=1:length(subs)

    %% 4.1. Reconstruction of source-level time-series

    % Read necessary files
    cd([outpath  subs{sub} '\' sess{sub} ])
    load dataclean1 % part 1 of the session
    load dataclean2 % part 2 of the session
    
    datacleans = {dataclean1, dataclean2};
    
    for dcleans =1:numel(datacleans)
        if dcleans == 1       % for the first part of the session
            load source_forward_10mm_1
            load source_inverse_10mm_1
        elseif dcleans == 2   % for the second part of the session
            load source_forward_10mm_2
            load source_inverse_10mm_2
        end

        % Reconstruct source-space data
        time         = datacleans{dcleans}.time{1};
        voxel_inside = find(source.inside==1);
        Nvox         = length(voxel_inside);
        datasource   = zeros(Nvox,length(time));
        for i = 1:Nvox
            disp([num2str(i) ' / ' num2str(length(voxel_inside))])
            datasource(i,:) = source.avg.filter{voxel_inside(i)} * datacleans{dcleans}.trial{1};
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

        if dcleans == 1
            % Use only the first 2,5 minutes recording of dataclean1 for all participants
            if time(end) > 140      % >2,5 min
                t1 = findbin(time,140);
                t2 = findbin(time,time(end));
                time(t1:t2) = [];
                datasource(:,t1:t2) = [];
            end
            % Use only 2,5 minutes recording of dataclean2 for all participants
        elseif dcleans == 2
            if time(end) > 290      % >2,5 min
                t1 = findbin(time,290);
                t2 = findbin(time,time(end));
                time(t1:t2) = [];
                datasource(:,t1:t2) = [];
            end
        end

        % Organize source-reconstructed data in a new fieldtrip structure
        datacleans{dcleans}.trial{1}   = single(datasource);
        datacleans{dcleans}.time{1}    = time;
        datacleans{dcleans}.label      = {};
        datacleans{dcleans}.sampleinfo = [1 size(datacleans{dcleans}.trial{1},2)];
        for i = 1:Nvox
            datacleans{dcleans}.label{i} = ['V' num2str(i)];
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
                datacleans{dcleans} = ft_redefinetrial(cfg,datacleans{dcleans});
            end
        end


        %% 4.2. Frequency analysis parameters: foi, toi, t_ftimwin

        f         = 0.55:0.05:3.55;        % originally up to 4.6 (now up to 34.8 Hz)
        foi       = exp(f);                % logarithmically spaced frequencies
        t_ftimwin = 5./foi;                % time-window adapted to each frequency (5 cycles)
        t_fstep   = 0.1;                   % sliding time-window (100 ms steps)

        bd   = t_ftimwin(1).*2;            % remove borders = 2*time-window
        bdpt = bd.*datacleans{dcleans}.fsample;

        Ntr  = 0;
        toi2 = {};
        ct   = 1;
        valid_tr = [];
        for tr = 1:length(datacleans{dcleans}.trial)
            toi = bd:t_fstep:datacleans{dcleans}.time{tr}(end)-bd;
            if length(toi) > 0
                toi2{ct} = toi;
                Ntr = Ntr+length(toi);
                valid_tr(ct) = tr;
                ct = ct+1;
            end
        end


        %% 4.3. Frequency analysis computation

        Nvox  = length(datacleans{dcleans}.label);
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
            freq           = ft_freqanalysis(cfg, datacleans{dcleans});

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
        Ntr_ok=length(trok) % all number of trials/spectra of the subj

        cd([outpath  subs{sub} '\' sess{sub} ])
        save(['freq_allvox_10mm_steps100ms_' num2str(dcleans) '.mat'], 'powsp', 'foi', 'Ntr_ok');

    end
end

