%% Fingerprinting of within-session group in Arana et al. (2025).

% This script is employed for the identification of subjects according to their natural frequencies fingerprints and conduct analysis to these results.
% It performs the correlations, the matches, accuracy of identification, the descriptive statistics, the bootsrapping to obtain de confidence intervals, and the permutation analysis for highest accuracy per chance cahnging participants' labels.
% The output is Corr_matrix (matrix with all correlation coefficients -Kendall's tau), Matches_01 (matrix with matches 1, and non-matches 0), accuracy_real, descriptive_stats, bootstrapping_analysis and permutation_analysis.
% It also performs the plots for the correlation matrix and matches (outputs map_matches_parula, map_matches_reedgreen).
% The input is, for each subject, the vector of natural frequencies per voxel for the first part of the session (or fisrt map) and the vector for the second part of the session (second map).

clear all
close all
clc

restoredefaultpath
addpath ('Z:\Toolbox\fieldtrip-20230118');
ft_defaults
addpath(genpath('Z:\Fingerprinting\scripts\Final'));


dpath = 'G:\Fingerprinting\Omega_data';
outpath = 'G:\Fingerprinting\Results\WS_group';

% Subjects and their single session
subs = {'sub-0001' 'sub-0002' 'sub-0003' 'sub-0004' 'sub-0005' 'sub-0006' 'sub-0007' 'sub-0008' 'sub-0009' 'sub-0011' 'sub-0012' 'sub-0014' 'sub-0015' 'sub-0016' 'sub-0018' 'sub-0019' 'sub-0020' 'sub-0021' 'sub-0022' 'sub-0023' 'sub-0024' 'sub-0025' 'sub-0026' 'sub-0027' 'sub-0028' 'sub-0029' 'sub-0030' 'sub-0031' 'sub-0032' 'sub-0033' 'sub-0034' 'sub-0035' 'sub-0037' 'sub-0039' 'sub-0040' 'sub-0041' 'sub-0042' 'sub-0044' 'sub-0045' 'sub-0046' 'sub-0047' 'sub-0048' 'sub-0049' 'sub-0050' 'sub-0051' 'sub-0052' 'sub-0055' 'sub-0056' 'sub-0057' 'sub-0058' 'sub-0059' 'sub-0060' 'sub-0061' 'sub-0062' 'sub-0063' 'sub-0064' 'sub-0065' 'sub-0067' 'sub-0068' 'sub-0069' 'sub-0070' 'sub-0071' 'sub-0072' 'sub-0073' 'sub-0074' 'sub-0075' 'sub-0076' 'sub-0077' 'sub-0078' 'sub-0079' 'sub-0080' 'sub-0084' 'sub-0085' 'sub-0087' 'sub-0088' 'sub-0089' 'sub-0090' 'sub-0091' 'sub-0092' 'sub-0094' 'sub-0095' 'sub-0096' 'sub-0097' 'sub-0098' 'sub-0099' 'sub-0101' 'sub-0102' 'sub-0103' 'sub-0104' 'sub-0105' 'sub-0106' 'sub-0134' 'sub-0145' 'sub-0146' 'sub-0148' 'sub-0149' 'sub-0150' 'sub-0151' 'sub-0152' 'sub-0154' 'sub-0155' 'sub-0156' 'sub-0157' 'sub-0158' 'sub-0159' 'sub-0160' 'sub-0161' 'sub-0165' 'sub-0166' 'sub-0167' 'sub-0168' 'sub-0169' 'sub-0170' 'sub-0171' 'sub-0175' 'sub-0176' 'sub-0177' 'sub-0179' 'sub-0181' 'sub-0184' 'sub-0185' 'sub-0195' 'sub-0197' 'sub-0200' 'sub-0207' 'sub-0208' 'sub-0210' 'sub-0212'};
sess = {'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0002' 'ses-0003' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0003' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001'};

% Natural frequencies from each half of the session
fnats1 = 'singlesub_fnat_steps200ms_part1';
fnats2 = 'singlesub_fnat_steps200ms_part2';

n_voxels = 1925; % number of voxels

% List all participants
AllParticipants = subs;
N_participants = length(AllParticipants); % total number of participants

% 7.1. Calculate correlations (output: Corr_matrix, Matches_01)
% 7.2. Calculate accuracy (output:accuracy_real)
% 7.3. Obtain plots of matches (output: figures map_matches_parula, map_matches_reedgreen)
% 7.4. Descriptive statistics (output: mean_self_corr, std_self_corr, std_self_corr, h, p_value, stats)
% 7.5. Bootstrapping (output: CI_lower, CI_lower)
% 7.6. Permutation analysis (output: max_correct_matches, accuracy_perm, p_value)

%% 7.1. Calculate correlations

% Initialize matrices
first_maps = cell(N_participants, 1);
second_maps = cell(N_participants, 1);
match_correlations = zeros(N_participants, 1);
match_participants = cell(N_participants, 1);

% Read part 1 of the session
for i = 1:N_participants
    cd([dpath '\' subs{i} '\' sess{i}]);
    load(fnats1);
    first_maps{i} = fnat.fnat';
end

% Read part 2 of the session
Corr_matrix = zeros(N_participants); % to save correlations
for i = 1:N_participants
    cd([dpath '\' subs{i} '\' sess{i}]);
    load(fnats2);
    second_maps{i} = fnat.fnat';

    for j = 1:N_participants
        % Correlation between maps
        Correlation_coefs = corr(first_maps{j}, second_maps{i}, 'Type', 'Kendall');
        Corr_matrix(i, j) = Correlation_coefs(1, 1); % save correlation
    end
end

% Identify matches
Matches = zeros(N_participants);
for sub = 1:size(Corr_matrix, 1)
    [corr_val, ind_f] = max(Corr_matrix(sub, :)); % higher correlation for each participant
    Matches(sub, ind_f) = corr_val;
    match_correlations(sub) = corr_val; % value of correlation
    match_participants{sub, 1} = AllParticipants{ind_f}; % id of participant
    match_participants{sub, 2} = ind_f; % number of participant
end
Matches_01=ceil(Matches); % matrix of matches (0 non-match, 1 match)

cd(outpath);
save Corr_matrix Corr_matrix;
save Matches_01 Matches_01

%% 7.2. Calculate accuracy
n_correct_matches = sum(diag(Matches_01)); % number of matches in diagonal (correct matches)
accuracy_real = (n_correct_matches / N_participants)*100;

cd(outpath);
save accuracy_real accuracy_real;

%% 7.3. Obtain plots of matches

% Matches parula
matches_parula = figure;
imagesc(Corr_matrix);

caxis([-0.8 0.8]); % adjust according to data (values of correlations)
cb = colorbar('eastoutside');

% Colorbar ticks
ticks = [-0.8 -0.4 0 0.4 0.8];
set(cb, 'Ticks', ticks);

colormap(parula); % plot values in parula

hold on;
[n_rows, n_cols] = size(Corr_matrix);


xlim([0.5, N_participants+0.5]);
xticks([0, 200]); % adjust according to data (N_participants)
ylim([0.5, N_participants+0.5]);
yticks([0, 200]); % adjust according to data (N_participants)
hold off;

cd(outpath);
saveas(matches_parula, 'map_matches_parula', 'png')

% Matches red-green
matches_red_green = figure;
plot_matches_redgreen(Matches); % pers function to map matches in red and green

cd(outpath);
saveas(matches_red_green, 'map_matches_reedgreen', 'png')

%% 7.4. Descriptive statistics

% Extract self-correlations (main diagonal)
self_correlations = diag(Corr_matrix);

% Initialize a vector to store the mean correlations with others for each subject
mean_non_self_correlations = zeros(N_participants, 1);

% Calculate the mean correlations with others for each subject
for i = 1:N_participants
    % Exclude the diagonal (self-correlation)
    non_self_values = Corr_matrix(i, :);
    non_self_values(i) = []; % Remove self-correlation
    % Calculate the mean of correlations with others
    mean_non_self_correlations(i) = mean(abs(non_self_values)); % abs to set positive the negative correlations
end

% Calculate descriptive statistics for self-correlations
mean_self_corr = mean(abs(self_correlations)); % abs to set positive the negative correlations 
std_self_corr = std(abs(self_correlations));

% Calculate descriptive statistics for correlations with others
mean_non_self_corr = mean(mean_non_self_correlations);
std_non_self_corr = std(mean_non_self_correlations);

% Display descriptive results
fprintf('Mean of self-correlations: %.4f\n', mean_self_corr);
fprintf('Standard deviation of self-correlations: %.4f\n', std_self_corr);
fprintf('Mean of non-self correlations: %.4f\n', std_non_self_corr);
fprintf('Standard deviation of non-self correlations: %.4f\n', std_non_self_corr);

% Perform paired t-test
[h, p_value, ci, stats] = ttest(self_correlations, mean_non_self_correlations);

% Display t-test results
fprintf('Paired t-test p-value: %.4f\n', p_value);
if h == 1
    fprintf('The self-correlations and mean non-self correlations are significantly different.\n');
else
    fprintf('There is insufficient evidence to say that the self-correlations and mean non-self correlations are significantly different.\n');
end

cd(outpath);
save descriptive_stats mean_self_corr std_self_corr std_self_corr h p_value stats

%% 7.5. Bootstrapping
rng('default')
n_bootstraps = 1000; % number of bootstraps
bootstrap_accuracies = zeros(n_bootstraps, 1); 

for b = 1:n_bootstraps
    % Resample participants with replacement
    resample_indices = randsample(N_participants, N_participants, true);
    resampled_Corr_matrix = Corr_matrix(resample_indices, resample_indices);
    
    % Identify correct matches
    Matches_boot = zeros(N_participants, 1);
    for sub = 1:N_participants
        [~, ind_f_boot] = max(resampled_Corr_matrix(sub, :)); % Higher correl after bootstrap
        if resample_indices(sub) == resample_indices(ind_f_boot)
            Matches_boot(sub) = 1; % If correct match
        end
    end
    
    % Calculate accuracy 
    bootstrap_accuracies(b) = (sum(Matches_boot) / N_participants) * 100;
end

CI_lower = prctile(bootstrap_accuracies, 2.5); %  95%
CI_upper = prctile(bootstrap_accuracies, 97.5); % 95%

fprintf('CI at 95%% es [%0.2f, %0.2f]\n', CI_lower, CI_lower);

cd(outpath);
save bootstrapping_analysis CI_lower CI_lower

%% 7.6. Permutation analysis
n_permutations = 1000; % number of perms

% Initialize
accuracy_perm = zeros(1, n_permutations);
n_correct_matches_perm = zeros(1, n_permutations);
rng('default')

for iter = 1:n_permutations
    % Permute labels of second map
    permuted_indices = randperm(N_participants);
    permuted_Corr_matrix = Corr_matrix(:, permuted_indices);

    % Calculate accuracy for perm
    Matches_perm = zeros(N_participants);
    for sub = 1:N_participants
        [corr_val, ind_f] = max(permuted_Corr_matrix(sub, :));
        Matches_perm(sub, ind_f) = corr_val;
    end

    n_correct_matches_perm = sum(diag(Matches_perm) > 0);
    accuracy_perm(iter) = n_correct_matches_perm / N_participants;
end


% Calculate p-value
max_correct_matches = max(n_correct_matches_perm);
p_value = sum(n_correct_matches_perm >= n_correct_matches) / n_permutations;
accuracy_perm = (max_correct_matches/N_participants)*100;

fprintf('Real accuracy: %.4f\n', accuracy_real);
fprintf('P-value of perm: %.4f\n', p_value);

message = sprintf('Across %d iterations, the highest success rate achieved was %d/%d, or  %.2f%%. Thus the p value associated with obtaining at least %d correct identifications is %.2f.', ...
    n_permutations, max_correct_matches, N_participants, accuracy_perm, n_correct_matches, p_value);

disp(message);

cd(outpath);
save permutation_analysis max_correct_matches accuracy_perm p_value
