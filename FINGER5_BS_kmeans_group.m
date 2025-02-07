%% K-means of between-session group in Arana et al. (2025).

% This script is employed for obtention of the power spectra centroids or clusters.
% Specifically, this was employed for k-means clustering of power spectra, of the between-session group but only from the fisrt session of each subject.
% 150 power spectra randomly selected of the first session per voxel and subject was introduced.
% The main output is the proportion of power spectra per cluster and the peak value of the centroid saved in kmeans_10mm_Nk25_150tr_step200ms.
% The input are the power spectra of the first session per subject (freq_allvox_10mm_steps100ms).
% Computation of K-means clustering requires near 32 Gb of RAM for the complete sample.

clear all
close all
clc

restoredefaultpath
addpath ('Z:\Toolbox\fieldtrip-20230118');
ft_defaults
addpath(genpath('Z:\Fingerprinting\scripts\Final'));


dpath = 'Z:\OMEGA\OMEGA_data\';
outpath = 'G:\Fingerprinting\kmeans\';


% Subjects and first session of between-session group
subs = {'sub-0001'	'sub-0002'	'sub-0006'	'sub-0008'	'sub-0011'	'sub-0016'	'sub-0019'	'sub-0020'	'sub-0022'	'sub-0023'	'sub-0025'	'sub-0030'	'sub-0032'	'sub-0035'	'sub-0039'	'sub-0040'	'sub-0041'	'sub-0042'	'sub-0044'	'sub-0046'	'sub-0048'	'sub-0049'	'sub-0050'	'sub-0051'	'sub-0106'	'sub-0150'	'sub-0200'};
sess = {'ses-0001'	'ses-0004'	'ses-0001'	'ses-0001'	'ses-0002'	'ses-0002'	'ses-0003'	'ses-0004'	'ses-0002'	'ses-0001'	'ses-0003'	'ses-0001'	'ses-0001'	'ses-0003'	'ses-0001'	'ses-0005'	'ses-0003'	'ses-0001'	'ses-0002'	'ses-0002'	'ses-0004'	'ses-0002'	'ses-0004'	'ses-0001'	'ses-0001'	'ses-0001'	'ses-0001'};


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
    load freq_allvox_10mm_steps100ms

    powsp = powsp(:,:,1:2:end); % select spectra each 200 ms steps


    ct = 1;
    valid_tr = [];

    for tr = 1:size(powsp,3)
        if sum(sum(isnan(powsp(:,:,tr)))) == 0
            valid_tr(ct) = tr;
            ct = ct+1;
        end
    end

    Ntr      = 150; % number of spectra per subject into kmeans

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
save kmeans_10mm_powsp_200ms_bs_group powsptot ksub kvox ktrial;


%% 5.2. K-means clustering
% Computation of K-means clustering
% Proportion of power spectra (out of 150) categorized in each cluster

% Read matrix with the power spectra of all subjects, time segments and voxels
Nk = 25;


% K-means clustering
[idx,C,sumd,D] = kmeans(powsptot,Nk,'Distance','cosine','Display','iter','Replicates',5,'MaxIter',200);


%% 5.3. Proportion of power spectra (out of 150) categorized in each cluster
% This computation is done for each subject and voxel (for each
% subject and voxel, the sum across clusters is equal to 1)

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
save kmeans_10mm_Nk25_150tr_step200ms_bs_group idx C sumd D Nvox propk -v7.3

