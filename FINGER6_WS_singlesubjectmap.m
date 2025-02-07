%% Single subject maps of within-session group in Arana et al. (2025).

% This script is employed for the single subject mapping of natural frequencies. 
% Specifically, this was employed for subjects within-session group that uses the both halves of one MEG session for each subject, separately. 
% The output is a natural frequencies vector of 1925 values for each half
% (singlesub_fnat_steps200ms_part1 and singlesub_fnat_steps200ms_part2) and both figures of each map.
% The inputs are freq_allvox_10mm_steps100ms_1, freq_allvox_10mm_steps100ms_2, kmeans_10mm_Nk25_150tr_step200ms_ws_group and sources for plotting.

clear all
close all
clc

restoredefaultpath
addpath ('Z:\Toolbox\fieldtrip-20230118');
ft_defaults
addpath ('Z:\Fingerprinting\scripts');

dpath = 'G:\Fingerprinting\Omega_data\';
figpath = 'G:\Fingerprinting\kmeans\Subject_maps_ws_group\' ;


% Subs, sessions and power sectra of both parts of the session of within-session group
subs = {'sub-0001' 'sub-0002' 'sub-0003' 'sub-0004' 'sub-0005' 'sub-0006' 'sub-0007' 'sub-0008' 'sub-0009' 'sub-0011' 'sub-0012' 'sub-0014' 'sub-0015' 'sub-0016' 'sub-0018' 'sub-0019' 'sub-0020' 'sub-0021' 'sub-0022' 'sub-0023' 'sub-0024' 'sub-0025' 'sub-0026' 'sub-0027' 'sub-0028' 'sub-0029' 'sub-0030' 'sub-0031' 'sub-0032' 'sub-0033' 'sub-0034' 'sub-0035' 'sub-0037' 'sub-0039' 'sub-0040' 'sub-0041' 'sub-0042' 'sub-0044' 'sub-0045' 'sub-0046' 'sub-0047' 'sub-0048' 'sub-0049' 'sub-0050' 'sub-0051' 'sub-0052' 'sub-0055' 'sub-0056' 'sub-0057' 'sub-0058' 'sub-0059' 'sub-0060' 'sub-0061' 'sub-0062' 'sub-0063' 'sub-0064' 'sub-0065' 'sub-0067' 'sub-0068' 'sub-0069' 'sub-0070' 'sub-0071' 'sub-0072' 'sub-0073' 'sub-0074' 'sub-0075' 'sub-0076' 'sub-0077' 'sub-0078' 'sub-0079' 'sub-0080' 'sub-0084' 'sub-0085' 'sub-0087' 'sub-0088' 'sub-0089' 'sub-0090' 'sub-0091' 'sub-0092' 'sub-0094' 'sub-0095' 'sub-0096' 'sub-0097' 'sub-0098' 'sub-0099' 'sub-0101' 'sub-0102' 'sub-0103' 'sub-0104' 'sub-0105' 'sub-0106' 'sub-0134' 'sub-0145' 'sub-0146' 'sub-0148' 'sub-0149' 'sub-0150' 'sub-0151' 'sub-0152' 'sub-0154' 'sub-0155' 'sub-0156' 'sub-0157' 'sub-0158' 'sub-0159' 'sub-0160' 'sub-0161' 'sub-0165' 'sub-0166' 'sub-0167' 'sub-0168' 'sub-0169' 'sub-0170' 'sub-0171' 'sub-0175' 'sub-0176' 'sub-0177' 'sub-0179' 'sub-0181' 'sub-0184' 'sub-0185' 'sub-0195' 'sub-0197' 'sub-0200' 'sub-0207' 'sub-0208' 'sub-0210' 'sub-0212'};
sess = {'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0002' 'ses-0003' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0003' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001'};
parts = {'freq_allvox_10mm_steps100ms_1.mat' 'freq_allvox_10mm_steps100ms_2.mat'};


cd('G:\Fingerprinting\kmeans\'); % load group Nk according to condition!!!
load kmeans_10mm_Nk25_150tr_step200ms_ws_group % kmeans for within-session group

cd('G:\Fingerprinting\Omega_data\sub-0001\ses-0001') % load sources of a subject as template to conduct the plotting
load source_forward_10mm_1
load source_inverse_10mm_1

voxel_inside = find(source.inside==1); % voxel inside cortical mask
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

    for part = 1:numel(parts)
        route = fullfile(dpath, subs{sub}, sess{sub}); % reconstruct route
        cd (route);
        load (parts{part});

        %% 6.1. Preparation of power spectra matrix for subsequent K-means clustering
        
        powsp = powsp(:,:,1:2:end); % to select all spectra each 200ms steps

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
       

        %% 6.2. Classification of power spectra into previous clusters of the whole group

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

        for vx=1:Nvox % smooth with neighbours
            vxneigh = connmat(vx,:)==1;
            propkzneig= propkz(:,vxneigh);
            propinterp = [];
            for i=1:size(propkzneig,2)
                propinterp(:,i) = interp(propkzneig(:,i),10);
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

        save_filename = sprintf('singlesub_fnat_steps200ms_part%d.mat', part);
        save(save_filename, 'fnat', 'propkz');

        %% 6.5. Plots of single-subject natural frequency maps

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
        
        cd(figpath); % folder to save figures
        print_filename = sprintf('%s_%s_singlesub_finger_part%d.tiff', subs{sub}, sess{sub}, part);
        print('-dtiff','-r300', print_filename);
        close
    end
end

