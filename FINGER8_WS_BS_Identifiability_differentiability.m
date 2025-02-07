%% This script should be run after FINGER7_WS_fingerprintng or FINGER7_BS_fingerprintng to calculate Identifiability and Differentiability
% The input needed is the matrix of correlations obtained in FINGER7 (Corr_matrix) and Matches_01. 
% You can use this script with Corr_matrix and Matches_01 of within-session group or with Corr_matrix and Matches_01 of between-session group. Load the corresponding Corr_matrix and Matches_01 accordingly.
% The outputs are graphic bars of Differentiability and Identifiability and the statistics.

clear all
close all
clc

restoredefaultpath
addpath ('Z:\Toolbox\fieldtrip-20230118');
ft_defaults
addpath(genpath('Z:\Fingerprinting\scripts\Final'));

dpath = 'G:\Fingerprinting\Results\WS_group'; % folder with variable Corr_matrix depending if WS or BS
cd(dpath); 
load Corr-matrix
load Matches_01

Iself_individual = diag(Corr_matrix);  % Self correlations
n = length(Iself_individual);     % Number of subjects

% Initialize
Iothers_mean = zeros(n, 1);
Iothers_std = zeros(n, 1);  
Differentiability = zeros(n, 1);

% Calculate Iothers_mean and Iothers_std for each subject
for i = 1:n
    % Exclude diagonal
    Iothers = Corr_matrix(i, :);
    Iothers(i) = [];  
    
    % Calculate Iothers
    Iothers_mean(i) = mean(Iothers);
    Iothers_std(i) = std(Iothers);  
end

% Calculate differentiability σij
for i = 1:n
    Differentiability(i) = (Iself_individual(i) - Iothers_mean(i)) / Iothers_std(i);
end

%% Figure

% Create the figure
difer = figure;
% Adjust the size of the figure (width, height, x position, y position)
set(difer, 'Position', [100, 100, 1200, 800]); % [x, y, width, height]

% Create non-stacked (overlaid) bar chart for identifiability
subplot(2, 1, 1);
% Ensure the base of the bar is non-negative
Iothers_mean_corrected = max(Iothers_mean, 0);  

% Create separate bars for the two data sets
hold on;
hBar1 = bar(Iself_individual, 'FaceColor', [0.678, 0.847, 0.902], 'BarWidth', 0.8); % Light blue bar
hBar2 = bar(Iothers_mean_corrected, 'FaceColor', [0.627, 0.627, 0.627], 'BarWidth', 0.8); % Gray bar

% Add transparency to both bars
hBar2.FaceAlpha = 0.6; % 70% transparency for the gray bar

set(gca, 'box', 'on');

% Axis labels and limits
% xticks(1:10:n);  % Labels every 10 subjects
ylim([0, 1]);  % Ensure Y-axis is between 0 and 1
yticks(0:0.2:1);  % Y-axis tick marks every 0.2

% Create a bar chart for Differentiability in the second panel
subplot(2, 1, 2);

% Create the bar chart and specify the color
bar(Differentiability, 'FaceColor', [0.678, 0.847, 0.902]); % Light blue for differentiability

% Axis labels and limits
% xticks('mode')
ylim([0,6]);
yticks(0:6);  

% Enhance visualization
grid on;
hold off;
cd(dpath);
saveas(difer, 'diff', 'png');

%% Stats in Identifiability and Differentiability

Iself_individual = diag(Corr_matrix);  % Autocorrelation (main diagonal)
n = length(Iself_individual);     % Number of subjects

% Initialize matrices for Iothers_mean, Iothers_std, Differentiability, and Identifiability
Iothers_mean = zeros(n, 1);
Iothers_std = zeros(n, 1);  % Individual standard deviation
Differentiability = zeros(n, 1);
Identifiability = zeros(n, 1);  % New variable for identifiability

% Initialize lists for matches (autocorrelations) and no matches (between different subjects)
matches_diff = [];
no_matches_diff = [];
matches_ident = [];
no_matches_ident = [];

% Calculate Iothers_mean and Iothers_std for each subject
for i = 1:n
    % Exclude the subject's own autocorrelation (diagonal)
    Iothers = Corr_matrix(i, :);
    Iothers(i) = [];  % Remove autocorrelation value (diagonal)
    
    % Calculate the mean and standard deviation of Iothers (correlations with other subjects)
    Iothers_mean(i) = mean(Iothers);
    Iothers_std(i) = std(Iothers);  % Individual standard deviation
end

% Calculate differentiability and identifiability for each subject
for i = 1:n
    Differentiability(i) = (Iself_individual(i) - Iothers_mean(i)) / Iothers_std(i);  % z-normalized differentiability
    Identifiability(i) = Iself_individual(i) - Iothers_mean(i);  % Identifiability (Iself - Iothers_mean)
    
    % Classify into matches (diagonal) and no matches (off-diagonal) using Matches_01
    if Matches_01(i, i) == 1  % If it's a match on the diagonal
        matches_diff = [matches_diff; i, Differentiability(i)];  % Store Differentiability for matches
        matches_ident = [matches_ident; i, Identifiability(i)];  % Store Identifiability for matches
    end
    
    % Classify no matches (off-diagonal)
    for j = 1:n
        if i ~= j && Matches_01(i, j) == 1  % If it's 1 off the diagonal
            no_matches_diff = [no_matches_diff; i, Differentiability(i)];  % Store Differentiability for no matches
            no_matches_ident = [no_matches_ident; i, Identifiability(i)];  % Store Identifiability for no matches
        end
    end
end


% Perform independent t-test for Identifiability (Matches vs No Matches)
fprintf('Performing independent t-test for Identifiability between matches and no matches.\n');
[h_ident, p_ident, ci_ident, stats_ident] = ttest2(matches_ident(:, 2), no_matches_ident(:, 2));

% Results for Identifiability
fprintf('Independent t-test p-value for Identifiability: %.4f\n', p_ident);
fprintf('t-statistic for Identifiability: %.4f\n', stats_ident.tstat);
fprintf('Degrees of freedom for Identifiability: %.0f\n', stats_ident.df);
fprintf('Confidence Interval for Identifiability: [%.4f, %.4f]\n', ci_ident(1), ci_ident(2));


% Perform independent t-test for Differentiability (Matches vs No Matches)
fprintf('Performing independent t-test for Differentiability between matches and no matches.\n');
[h_diff, p_diff, ci_diff, stats_diff] = ttest2(matches_diff(:, 2), no_matches_diff(:, 2));

% Results for Differentiability
fprintf('Independent t-test p-value for Differentiability: %.4f\n', p_diff);
fprintf('t-statistic for Differentiability: %.4f\n', stats_diff.tstat);
fprintf('Degrees of freedom for Differentiability: %.0f\n', stats_diff.df);
fprintf('Confidence Interval for Differentiability: [%.4f, %.4f]\n', ci_diff(1), ci_diff(2));


% Calculate the mean and standard deviation of Iself_individual, Iothers_mean, Differentiability, and Identifiability
mean_Iself = mean(Iself_individual);          % Mean of Iself (autocorrelations)
std_Iself = std(Iself_individual);            % Standard deviation of Iself

mean_Iothers = mean(Iothers_mean);            % Mean of Iothers (correlations with other individuals)
std_Iothers = std(Iothers_mean);              % Standard deviation of Iothers

mean_Iothers_std = mean(Iothers_std);         % Mean of the standard deviations of Iothers
std_Iothers_std = std(Iothers_std);           % Standard deviation of the standard deviations of Iothers

mean_Identifiability = mean(Identifiability);  % Mean of Identifiability
std_Identifiability = std(Identifiability);    % Standard deviation of Identifiability

mean_Differentiability = mean(Differentiability);  % Mean of Differentiability
std_Differentiability = std(Differentiability);    % Standard deviation of Differentiability


% Calculate the mean and standard deviation of Differentiability and Identifiability in matches and non_matches
mean_non_matches_dif = mean(no_matches_diff(:, 2));
std_non_matches_dif = std(no_matches_diff(:, 2));
mean_matches_dif = mean(matches_diff(:, 2));
std_matches_dif = std(matches_diff(:, 2));

mean_non_matches_ident = mean(no_matches_ident(:, 2));
std_non_matches_ident = std(no_matches_ident(:, 2));
mean_matches_ident = mean(matches_ident(:, 2));
std_matches_ident = std(matches_ident(:, 2));

