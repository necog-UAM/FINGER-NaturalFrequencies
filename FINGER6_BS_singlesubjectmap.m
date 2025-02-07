%% Single subject maps of between-session group in Arana et al. (2025).

% This script is employed for the single subject mapping of natural frequencies.
% Specifically, this was employed for subjects between-session group that uses the both sessions for each subject, separately.
% The output is a natural frequencies vector of 1925 for each session (singlesub_fnat_steps200m) ann a figure of the map of each session.
% The inputs are freq_allvox_10mm_steps100ms, kmeans_10mm_Nk25_150tr_step200ms_bs_group and sources for plotting.

clear all
close all
clc

restoredefaultpath
addpath ('Z:\Toolbox\fieldtrip-20230118');
ft_defaults
addpath(genpath('Z:\Fingerprinting\scripts\Final'));


dpath = 'Z:\OMEGA\OMEGA_data\';
figpath = 'G:\Fingerprinting\kmeans\Subject_maps_bs_group\' ;


% Subs and sessions of between-session group
subs = {'sub-0001'	'sub-0001'	'sub-0002'	'sub-0002'	'sub-0006'	'sub-0006'  'sub-0008'	'sub-0008'	'sub-0011'	'sub-0011'  'sub-0016'	'sub-0016'	'sub-0019'	'sub-0019'	'sub-0020'	'sub-0020'	'sub-0022'	'sub-0022'	'sub-0023'	'sub-0023'	'sub-0025'	'sub-0025'	'sub-0030'	'sub-0030'	'sub-0032'	'sub-0032'	'sub-0035'	'sub-0035'	'sub-0039'	'sub-0039'	'sub-0040'	'sub-0040'	'sub-0041'	'sub-0041'	'sub-0042'	'sub-0042'	'sub-0044'	'sub-0044'	'sub-0046'	'sub-0046'	'sub-0048'	'sub-0048'	'sub-0049'	'sub-0049'	'sub-0050'	'sub-0050'	'sub-0051'	'sub-0051'	'sub-0106'	'sub-0106'	'sub-0150'	'sub-0150'	'sub-0200'	'sub-0200'};
sess = {'ses-0001'	'ses-0003'	'ses-0004'	'ses-0002'	'ses-0001'	'ses-0003'  'ses-0001'	'ses-0002'	'ses-0002'	'ses-0001'  'ses-0002'	'ses-0001'	'ses-0003'	'ses-0002'	'ses-0004'	'ses-0003'	'ses-0002'	'ses-0001'	'ses-0001'	'ses-0002'	'ses-0003'	'ses-0002'	'ses-0001'	'ses-0002'	'ses-0001'	'ses-0002'	'ses-0003'	'ses-0002'	'ses-0001'	'ses-0002'	'ses-0005'	'ses-0004'	'ses-0003'	'ses-0001'	'ses-0001'	'ses-0002'	'ses-0002'	'ses-0003'	'ses-0002'	'ses-0001'	'ses-0004'	'ses-0002'	'ses-0002'	'ses-0001'	'ses-0004'	'ses-0002'	'ses-0001'	'ses-0002'	'ses-0001'	'ses-0002'	'ses-0001'	'ses-0002'	'ses-0001'	'ses-0002'};

cd('G:\Fingerprinting\kmeans\'); % load group Nk according to condition!!!
load kmeans_10mm_Nk25_150tr_step200ms_bs_group % kmeans for between-session group

cd('Z:\OMEGA\OMEGA_data\sub-0001\ses-0003'); % load sources of a subject as template to conduct the plotting
load source_forward_10mm
load source_inverse_10mm

voxel_inside = find(source.inside==1);  % voxel inside cortical mask
Nvox= 1925;
Nk = 25;


% 6.1. Preparation of power spectra matrix for subsequent K-means clustering
% 6.2. Classification of power spectra into previous clusters of the whole group
% 6.3. Proportion of power spectra categorized in each cluster
% 6.4. Compute natural frequency of each voxel (output: fnat and propkz)
% 6.5. Plots of single-subject natural frequency maps (output: figures)

for sub=1:numel(subs)

    route = fullfile(dpath, subs{sub}, sess{sub}); % reconstruct route
    cd (route);


    load freq_allvox_10mm_steps100ms

    powsp = powsp(:,:,1:2:end); % to select power spectra each 200ms steps

    %% 6.1. Preparation of power spectra matrix for subsequent K-means clustering


    rng('shuffle')

    ct = 1;
    powspsub  = single(zeros(size(powsp,1)*size(powsp,3),size(powsp,2)));
    kvox   = [];
    ktrial = [];
    for i = 1:size(powsp,1)
        for tr = 1:size(powsp,3)
            powspsub(ct,:) = powsp(i,:,tr);
            kvox(ct)    = i;
            ktrial(ct)  = tr;
            ct = ct+1;
        end
    end

    bl       = sum(powspsub,2);                                % save kmeans_10mm_baseline bl
    powspsub = powspsub./repmat(bl,[1 size(powspsub,2)]);      % compute relative power to correct the center of the head bias


    %% 6.2. Classification of power spectra into clusters of the whole group

    % K-means clustering
    [D,I] = pdist2(C,powspsub, 'cosine','Smallest',1);

    %% 6.3. Proportion of power spectra categorized in each cluster

    propk = NaN(Nk,Nvox);
    for k = 1:Nk
        disp(['Cluster ' num2str(k) '/' num2str(Nk)])
        ctk = find(I==k);
        ctkvox = kvox(ctk)';
        ctkvoxfilt = ctkvox==[1:Nvox];
        propk(k,:)=sum(ctkvoxfilt)./Ntr;
    end

    %% 6.4. Compute natural frequency of each voxel

    f   = 0.55:0.05:3.55;        % until 34.8 Hz
    foi = exp(f);
    ff  = [];
    for k = 1:Nk
        sp = C(k,:);
        [pks,locs] = findpeaks(sp);
        fx  = round(foi(locs(find(pks == max(pks)))),1);
        ff(k) = fx;
    end

    [ffsort,idf]=sort(ff);
    ffsort = unique(ffsort);

    ff2=ff;
    ff = NaN(Nk,2);
    badk = [];
    %     figure
    for k = 1:Nk
        sp = C(idf(k),:);
        %         plot(foi,sp)
        %         title(num2str(ff2(idf(k))))
        %         pause
        [pks,locs] = findpeaks(sp,'MinPeakHeight',0.1,'MinPeakProminence',0.02);
        % if it does not find any peak, go to the next centroid
        if length(pks) == 1                                % unimodal spectrum
            ff(idf(k),1) = round(foi(locs(1)),1);
            ff(idf(k),2) = round(foi(locs(1)),1);      % if unimodal spectrum, repeat 1st peak
        elseif length(pks) == 2                            % bimodal spectrum
            ff(idf(k),1) = round(foi(locs(1)),1);
            ff(idf(k),2) = round(foi(locs(2)),1);
        elseif length(pks) > 2                         % in case of a third residual peak
            [~,id] = sort(pks,'descend');
            locs = locs(id(1:2));
            ff(idf(k),1) = round(foi(locs(1)),1);
            ff(idf(k),2) = round(foi(locs(2)),1);
        elseif length(pks) == 0
            badk = [badk idf(k)];
        end
    end

    propk2 = zeros(length(ffsort),Nvox);

    for i=1:length(ffsort)
        [r,c] = find(ff==ffsort(i));
        r = unique(r);
        for j=1:length(r)
            if ff(r(j),1) ~= ff(r(j),2)
                propk2(i,:) = propk2(i,:) + 1/2.*sum(propk(r(j),:),1);    % bimodal
            elseif ff(r(j),1) == ff(r(j),2)
                propk2(i,:) = propk2(i,:) + sum(propk(r(j),:),1);    % unimodal
            end
        end
    end


    emptycol = find(sum(propk2,2)==0);
    propk2(emptycol,:) = [];
    ffsort(emptycol) = [];

    propkz = (propk2-mean(propk2,2))./std(propk2,0,2);        % z-normalize

    [dim,xx,yy,zz,connmat, dtempl] = Omega_neighbors_ly(source); % pers function

    propkz2sm = NaN(size(propkz));
    fnatsm=[];
    fnatsm_T=[];
    fnatsm_p=[];
    f2 = interp(ffsort,10);

    for vx=1:Nvox  % smooth with neighbours
        vxneigh = connmat(vx,:)==1;
        propkzneig= propkz(:,vxneigh);
        propinterp = [];
        for i=1:size(propkzneig,2)
            propinterp(:,i) = interp(propkzneig(:,i),10);
            % propinterp(:,i) = interp1(propkzneig(:,i),f2,'pchip');
        end
        [h,p,ci,stats] = ttest(propinterp');

        [tmax,idmax] = max(stats.tstat);
        fnatsm(vx) = f2(idmax);
        fnatsm_T(vx) = tmax;
        fnatsm_p(vx) = p(idmax);

        tval = stats.tstat;
        tval(p>.05) = NaN;
        tval(tval < 0) = NaN;
        [pks,locs] = findpeaks(tval);

        fnat.pks{vx}=f2(locs);
        fnat.w{vx}=pks*100./sum(pks);

    end

    fnatsm2 = fnatsm;
    fnatsm2(fnatsm_p>.05) = NaN;

    fnat.fnat = fnatsm;
    fnat.fnatsig = fnatsm2;
    fnat.tval = fnatsm_T;
    fnat.pval = fnatsm_p;
    cd (route);
    save singlesub_fnat_steps200ms fnat propkz

    %% 6.5.  Plots of single-subject natural frequency maps

    source2 = source;
    source2.avg.pow(voxel_inside) = log(fnatsm2);
    source2.avg.mom = cell(length(source2.inside),1);

    cfg=[];
    cfg.parameter  = 'pow';
    cfg.downsample = 2;
    cfg.interpmethod = 'nearest';
    source_interp = ft_sourceinterpolate (cfg, source2, source_forward.mri);

    figure('WindowState','maximized','Color',[1 1 1]);
    % figure

    cfg               = [];
    cfg.figure        = 'gca';
    cfg.method        = 'surface';
    cfg.funparameter  = 'pow';
    cfg.maskparameter = cfg.funparameter;
    cfg.funcolorlim   = [0.7 3.4];
    cfg.funcolormap   = 'jet_omega_mod';
    cfg.projmethod    = 'nearest';
    cfg.opacity       = 0.8;
    cfg.camlight      = 'no';
    cfg.colorbar      = 'no';
    cfg.surffile     = 'surface_pial_left.mat';
    cfg.surfinflated  = 'surface_inflated_left_caret_white.mat';
    subplot(2,2,1), ft_sourceplot(cfg,source_interp), view([-90 0]), camlight('left')
    subplot(2,2,3), ft_sourceplot(cfg,source_interp), view([90 0]),  camlight('left')

    cfg.surffile     = 'surface_pial_right.mat';
    cfg.surfinflated  = 'surface_inflated_right_caret_white.mat';
    subplot(2,2,2), ft_sourceplot(cfg,source_interp), view([90 0]),  camlight('right')
    subplot(2,2,4), ft_sourceplot(cfg,source_interp), view([-90 0]), camlight('right')

    cd(figpath); % folder to save figs
    print('-dtiff','-r300',[subs{sub} '_' sess{sub} '_singlesub_finger.tiff']);
    close
end

