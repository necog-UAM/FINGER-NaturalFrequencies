%% K-means of within-session group in Arana et al. (2025).

% This script is employed for obtention of the power spectra centroids or clusters.
% Specifically, this was employed for k-means clustering of power spectra, of the within-session group but only from the fisrt part of the session of each subject. 
% 150 power spectra randomly selected of the first part of the session per voxel and subject was introduced.
% The main output is the proportion of power spectra per cluster and the peak value of the centroid saved in kmeans_10mm_Nk25_150tr_step200ms.
% The input are the power spectra of the first part of the session per subject (freq_allvox_10mm_steps100ms_1).
% Computation of K-means clustering requires near 32 Gb of RAM for the complete sample.

clear all
close all
clc

restoredefaultpath
addpath ('Z:\Toolbox\fieldtrip-20230118');
ft_defaults
addpath(genpath('Z:\Fingerprinting\scripts\Final'));

dpath = ' G:\Fingerprinting\Omega_data'; 
outpath = 'G:\Fingerprinting\kmeans\';


% Subs and sess of within-session group.
subs = {'sub-0001' 'sub-0002' 'sub-0003' 'sub-0004' 'sub-0005' 'sub-0006' 'sub-0007' 'sub-0008' 'sub-0009' 'sub-0011' 'sub-0012' 'sub-0014' 'sub-0015' 'sub-0016' 'sub-0018' 'sub-0019' 'sub-0020' 'sub-0021' 'sub-0022' 'sub-0023' 'sub-0024' 'sub-0025' 'sub-0026' 'sub-0027' 'sub-0028' 'sub-0029' 'sub-0030' 'sub-0031' 'sub-0032' 'sub-0033' 'sub-0034' 'sub-0035' 'sub-0037' 'sub-0039' 'sub-0040' 'sub-0041' 'sub-0042' 'sub-0044' 'sub-0045' 'sub-0046' 'sub-0047' 'sub-0048' 'sub-0049' 'sub-0050' 'sub-0051' 'sub-0052' 'sub-0055' 'sub-0056' 'sub-0057' 'sub-0058' 'sub-0059' 'sub-0060' 'sub-0061' 'sub-0062' 'sub-0063' 'sub-0064' 'sub-0065' 'sub-0067' 'sub-0068' 'sub-0069' 'sub-0070' 'sub-0071' 'sub-0072' 'sub-0073' 'sub-0074' 'sub-0075' 'sub-0076' 'sub-0077' 'sub-0078' 'sub-0079' 'sub-0080' 'sub-0084' 'sub-0085' 'sub-0087' 'sub-0088' 'sub-0089' 'sub-0090' 'sub-0091' 'sub-0092' 'sub-0094' 'sub-0095' 'sub-0096' 'sub-0097' 'sub-0098' 'sub-0099' 'sub-0101' 'sub-0102' 'sub-0103' 'sub-0104' 'sub-0105' 'sub-0106' 'sub-0134' 'sub-0145' 'sub-0146' 'sub-0148' 'sub-0149' 'sub-0150' 'sub-0151' 'sub-0152' 'sub-0154' 'sub-0155' 'sub-0156' 'sub-0157' 'sub-0158' 'sub-0159' 'sub-0160' 'sub-0161' 'sub-0165' 'sub-0166' 'sub-0167' 'sub-0168' 'sub-0169' 'sub-0170' 'sub-0171' 'sub-0175' 'sub-0176' 'sub-0177' 'sub-0179' 'sub-0181' 'sub-0184' 'sub-0185' 'sub-0195' 'sub-0197' 'sub-0200' 'sub-0207' 'sub-0208' 'sub-0210' 'sub-0212'};
sess = {'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0002' 'ses-0003' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0003' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001'};


% 5.1. Preparation of power spectra matrix for subsequent K-means clustering (output: powsptot ksub kvox ktrial)
% 5.2. K-means clustering (output: idx C sumd D)
% 5.3. Proportion of power spectra categorized in each cluster (output: propk)


%% 5.1. Preparation of power spectra matrix for subsequent K-means clustering

powsptot = [];
ksub     = [];
kvox     = [];
ktrial   = [];
rng('default')
rng('shuffle')


for sub=1:length(subs)
    route = fullfile([dpath, subs{sub} '\' sess{sub}]); % reconstruct route
    cd(route)
    load freq_allvox_10mm_steps100ms_1  % first part of the session

    powsp = powsp(:,:,1:2:end);  % select spectra each 200 ms steps


    ct = 1;
    valid_tr = [];

    for tr = 1:size(powsp,3)
        if sum(sum(isnan(powsp(:,:,tr)))) == 0
            valid_tr(ct) = tr;
            ct = ct+1;
        end
    end

    Ntr      = 150;  % number of spectra selected per subject into kmeans

    rndtr   = randperm(length(valid_tr),Ntr); % caution if Ntr > length(valid<-tr)
    select_tr = valid_tr(rndtr);
    powsp2  = powsp(:,:,select_tr);
    ct = 1;
    powsp3  = single(zeros(size(powsp2,1)*size(powsp2,3),size(powsp2,2)));
    ksub2   = [];
    kvox2   = [];
    ktrial2 = [];
    for i = 1:size(powsp2,1)
        for tr = 1:size(powsp2,3)
            powsp3(ct,:) = powsp2(i,:,tr);
            ksub2(ct)    = sub;
            kvox2(ct)    = i;
            ktrial2(ct)  = select_tr(tr);
            ct = ct+1;
        end
    end

    powsptot = [powsptot; powsp3];         % concatenation of power spectra
    ksub     = [ksub ksub2];               % keep track of subject, voxel and trial of each power spectrum
    kvox     = [kvox kvox2];
    ktrial   = [ktrial ktrial2];
end

clear freq powsp

bl       = sum(powsptot,2);                                % save kmeans_10mm_baseline bl
powsptot = powsptot./repmat(bl,[1 size(powsptot,2)]);      % compute relative power to correct the center of the head bias


cd(outpath);
save kmeans_10mm_powsp_200ms_ws_group powsptot ksub kvox ktrial;


%% 5.2. K-means clustering
% Computation of K-means clustering
% Proportion of power spectra (out of 150) categorized in each cluster

% Read matrix with the power spectra of all subjects, time segments and voxels
Nk = 25;


% K-means clustering
[idx,C,sumd,D] = kmeans(powsptot,Nk,'Distance','cosine','Display','iter','Replicates',5,'MaxIter',200);


%% 5.3. Proportion of power spectra (out of 150) categorized in each cluster
% This computation is done for each subject and voxel (for each subject and voxel, the sum across clusters is equal to 1)

Nvox = 1925;
Ntr = 150;
Nsub=numel(subs);

propk = NaN(Nk,Nvox,Nsub);

for k = 1:Nk
    disp(['Cluster ' num2str(k) '/' num2str(Nk)])
    ctk = find(idx==k);
    ctksub = ksub(ctk)';
    ctkvox = kvox(ctk)';
    ctksubfilt = ctksub==[1:Nsub];
    ctkvoxfilt = ctkvox==[1:Nvox];
    for s=1:Nsub
        propk(k,:,s)=sum(ctksubfilt(:,s) & ctkvoxfilt)./Ntr;
    end
end


cd(outpath);
save kmeans_10mm_Nk25_150tr_step200ms_ws_group idx C sumd D Nvox propk -v7.3 

