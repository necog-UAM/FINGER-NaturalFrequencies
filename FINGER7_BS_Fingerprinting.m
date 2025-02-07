%% Fingerprinting of between-session group in Arana et al. (2025).

% This script is employed for the identification of subjects according to their natural frequencies fingerprints and conduct analysis to these results.
% It performs the correlations, the matches, accuracy of identification, the descriptive statistics, the bootsrapping to obtain de confidence intervals, and the permutation analysis for highest accuracy per chance cahnging participants' labels.
% The output is Corr_matrix (matrix with all correlation coefficients -Kendall's tau), Matches_01 (matrix with matches 1, and non-matches 0), accuracy_real, descriptive_stats, bootstrapping_analysis and permutation_analysis.
% It also performs the plots for the correlation matrix and matches (outputs map_matches_parula, map_matches_reedgreen).
% Finaly, the script conduct the analysis of days elapsed between sessions, specific for this group. The ouput is Days_elapsed_vs_self_correlation_Pearson and figure of correlations, Welch_matches_vs_days and violin figure with days elapsed for matches and for non-matches.
% The input is, for each subject, the vector of natural frequencies per voxel for the first session (or fisrt map) and the vector for the second session (second map).

clear all
close all
clc

restoredefaultpath
addpath ('Z:\Toolbox\fieldtrip-20230118');
ft_defaults
addpath(genpath('Z:\Fingerprinting\scripts\Final'));


dpath = 'G:\Fingerprinting\Omega_data';
outpath = 'G:\Fingerprinting\Results\BS_group';

% Subjects and their first and second sessions
subs = {'sub-0001'	'sub-0002'	'sub-0006'	'sub-0008'	'sub-0011'	'sub-0016'	'sub-0019'	'sub-0020'	'sub-0022'	'sub-0023'	'sub-0025'	'sub-0030'	'sub-0032'	'sub-0035'	'sub-0039'	'sub-0040'	'sub-0041'	'sub-0042'	'sub-0044'	'sub-0046'	'sub-0048'	'sub-0049'	'sub-0050'	'sub-0051'	'sub-0106'	'sub-0150'	'sub-0200'};
sess1 = {'ses-0001' 'ses-0004' 'ses-0001'  'ses-0001'	'ses-0002'	'ses-0002'	'ses-0003'	'ses-0004'	'ses-0002'	'ses-0001'	'ses-0003'	'ses-0001'	'ses-0001'	'ses-0003'	'ses-0001'	'ses-0005'	'ses-0003'	'ses-0001'  'ses-0002'	'ses-0002'	'ses-0004'	'ses-0002'	'ses-0004'	'ses-0001'  'ses-0001'	'ses-0001'	'ses-0001'};
sess2 = {'ses-0003' 'ses-0002' 'ses-0003'	'ses-0002'	'ses-0001' 	'ses-0001'	'ses-0002'	'ses-0003'	'ses-0001'	'ses-0002'	'ses-0002'	'ses-0002'	'ses-0002'	'ses-0002'	'ses-0002'	'ses-0004'	'ses-0001'	'ses-0002'	'ses-0003'	'ses-0001'	'ses-0002'	'ses-0001'	'ses-0002'	'ses-0002'	'ses-0002'	'ses-0002'	'ses-0002'};

% Natural frequencies saved on each corresponding session folder
fnats = 'singlesub_fnat_steps200ms.mat';

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
% 7.7. Days between sessions analysis (output: Days_elapsed_vs_self_correlation_Pearson -r, p; Welch_matches_vs_days -h, p, stats, mean_days_match, std_days_match, mean_days_no_match, std_days_no_match; 
% and figures Correl_days_vs_self_values, Violin_days_matches)

%% 7.1. Calculate correlations

% Initialize matrices
first_maps = cell(N_participants, 1);
second_maps = cell(N_participants, 1);
match_correlations = zeros(N_participants, 1);
match_participants = cell(N_participants, 1);

% Read session 1
for i = 1:N_participants
    cd([dpath '\' subs{i} '\' sess1{i}]);
    load(fnats);
    first_maps{i} = fnat.fnat';
end

% Read session 2
Corr_matrix = zeros(N_participants); % to save correlations
for i = 1:N_participants
    cd([dpath '\' subs{i} '\' sess2{i}]);
    load(fnats);
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
xticks([0, 30]); % adjust according to data (N_participants)
ylim([0.5, N_participants+0.5]);
yticks([0, 30]); % adjust according to data (N_participants)
hold off;

cd(outpath);
saveas(matches_parula, 'map_matches_parula', 'png')

% Matches red-green
matches_red_green = figure;
plot_matches_redgreen(Matches); % pers function to map matches in red and green

cd(outpath);
saveas(matches_red_green, 'map_matches_reedgreen', 'png')

%% 7.3. Descriptive statistics

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

%% 7.4. Bootstrapping
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

%% 7.5. Permutation analysis
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

%% 7.7. Days between sessions analysis 
% To examine if self_correlation values or number of correct matches are related to days elapsed between sessions of each subject.

% Number of the participant, code and days elapsed between sessions of the participant
subject_days = {
    1, 'sub-0001', 1696;
    2, 'sub-0002', 7;
    3, 'sub-0006', 22;
    4, 'sub-0008', 1029;
    5, 'sub-0011', 263;
    6, 'sub-0016', 799;
    7, 'sub-0019', 148;
    8, 'sub-0020', 305;
    9, 'sub-0022', 275;
    10, 'sub-0023', 96;
    11, 'sub-0025', 127;
    12, 'sub-0030', 30;
    13, 'sub-0032', 6;
    14, 'sub-0035', 61;
    15, 'sub-0039', 0;
    16, 'sub-0040', 550;
    17, 'sub-0041', 258;
    18, 'sub-0042', 94;
    19, 'sub-0044', 817;
    20, 'sub-0046', 13;
    21, 'sub-0048', 115;
    22, 'sub-0049', 8;
    23, 'sub-0050', 29;
    24, 'sub-0051', 852;
    25, 'sub-0106', 373;
    26, 'sub-0150', 103;
    27, 'sub-0200', 1
    };

self_correlations = diag(Corr_matrix);

% 7.7.1 Analysis to examine correlation (Pearson) of days elapsed between sessions and self-correlation values: 

% First, the plot of self_correlation and days between sessions
colors = zeros(N_participants, 3); % RGB initialize

% Calculate coincidences
for i = 1:N_participants
    if Corr_matrix(i, i) == max(Corr_matrix(i, :))
        colors(i, :) = [0.2, 0.7, 0.3]; % Green [R, G, B] correct match
    else
        colors(i, :) = [0.9, 0.3, 0.2]; % Red [R, G, B] incorrect match
    end
end

% Extract days
days_between_sessions = zeros(N_participants, 1);
for i = 1:N_participants
    sub_index = find(strcmp(subject_days(:, 2), AllParticipants{i}));
    days_between_sessions(i) = subject_days{sub_index, 3};
end

% Graphic
days_correl = figure;
scatter(days_between_sessions, self_correlations, 50, colors, 'filled');
xlabel('Days Between Sessions');
ylabel('Self-Correlation (R)');
xlim([-60 1750]);
xticks([0 200 400 600 800 1000 1200 1400 1600 1800]);
ylim([0 0.8]);
grid off;

set(gcf, 'Position', [100, 100, 500, 500]);
set(gca, 'box', 'on');

% Linear
lm = fitlm(days_between_sessions, self_correlations, 'RobustOpts', 'on');

% Tendency line
x_fit = linspace(min(days_between_sessions), max(days_between_sessions), 100)';
y_fit = predict(lm, x_fit);

hold on;
plot(x_fit, y_fit, 'k-', 'LineWidth', 0.3);
hold off;

cd(outpath);
saveas(days_correl, 'Correl_days_vs_self_values', 'png')

% Now, calculate pearson coefficient for correlation between self_correlation values and days elapsed between sessions
[r, p] = corr(days_between_sessions, self_correlations, 'Type', 'Pearson');

fprintf('Pearson R-days: %.4f\n', r);
fprintf( 'p: %.4f\n', p);

cd(outpath);
save Days_elapsed_vs_self_correlation_Pearson r p


% 7.7.2 Examine the number of correct matches and non-correct macthes vs days between sessions:

% List of subjects that do NOT match with themselves (non-correct matches)
no_match_subjects = [2, 6, 10, 18, 22, 27];

% Initialize vectors for days between sessions for both groups
days_match = [];
days_no_match = [];

% Loop through all subjects and classify based on whether they match or not
for i = 1:N_participants
    if ismember(i, no_match_subjects)
        % If the subject is in the no-match list
        days_no_match = [days_no_match; subject_days{i, 3}];
    else
        % If the subject matches with themselves
        days_match = [days_match; subject_days{i, 3}];
    end
end

% Calculate the mean and standard deviation for correct matches and non-correct matches and days between sessions
mean_days_match = mean(days_match);
std_days_match = std(days_match);

mean_days_no_match = mean(days_no_match);
std_days_no_match = std(days_no_match);

% Display results
fprintf('Mean days between sessions for matching subjects: %.2f\n', mean_days_match);
fprintf('Standard deviation of days for matching subjects: %.2f\n', std_days_match);
fprintf('Mean days between sessions for non-matching subjects: %.2f\n', mean_days_no_match);
fprintf('Standard deviation of days for non-matching subjects: %.2f\n', std_days_no_match);


% Perform Welch's t-test
[h, p, ci, stats] = ttest2(days_match, days_no_match, 'Vartype', 'unequal');

% Display Welch's t-test results
fprintf('Results of Welch’s t-test:\n');
fprintf('t-value: %.4f\n', stats.tstat);
fprintf('Degrees of freedom: %.4f\n', stats.df);
fprintf('p-value: %.4f\n', p);

cd(outpath);
save Welch_matches_vs_days h p stats mean_days_match std_days_match mean_days_no_match std_days_no_match

% Now, violin plot of results:

% Data for days between sessions
Y = {days_match, days_no_match}; % Organize data into a cell array

% Create the violin plot
violin_days = figure;
[h, L, MX, MED] = violin(Y, 'xlabel', {'Correct match', 'Non-correct match'}, ...
    'facecolor', [0.2, 0.7, 0.3; 0.9, 0.3, 0.2], ... % Custom colors
    'edgecolor', 'k', ... % Border color for violins
    'facealpha', 0.3, ... % Transparency
    'plotlegend', 0, ...
    'mc', 'k', ... % Mean color
    'medc', 'b'); % Median color

% Add points inside the violins
hold on; % Keep current plot

% Jitter settings for points
jitter_amount = 0.1; % Amount of shift for points

% Draw points for 'Correct match'
scatter(ones(size(days_match)) + (rand(size(days_match)) - 0.5) * jitter_amount, ...
    days_match, 20, [0.2, 0.7, 0.3], 'filled', 'MarkerEdgeColor', 'none'); % Point color and style

% Draw points for 'Non-correct match'
scatter(2 * ones(size(days_no_match)) + (rand(size(days_no_match)) - 0.5) * jitter_amount, ...
    days_no_match, 20, [0.9, 0.3, 0.2], 'filled', 'MarkerEdgeColor', 'none'); % Point color and style

set(violin_days, 'Position', [100, 100, 500, 500]); % [x, y, width, height]
ylim([-400, 2100]);
yticks(0:200:1800); % Y-axis tick points every 2

% Add labels
ylabel('Days between sessions');
title('Comparison of Days between Sessions: Correct match vs Non-correct Match');

% Save plot as PNG file
cd(outpath);
saveas(violin_days, 'Violin_days_matches', 'png');

hold off; 
