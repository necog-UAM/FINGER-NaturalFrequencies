cd('Z:\OMEGA\OMEGA_data\sub-0001\ses-0001');
load source_forward_10mm
load source_inverse_10mm

% load first_maps
% load second_maps

N_participants = numel(first_maps);
n_voxels = 1925;

% 1. Compute baseline
Iself_real = zeros(N_participants,1);
Identifiability_real = zeros(N_participants,1);
Differentiability_real = zeros(N_participants,1);

for i = 1:N_participants
    self_corr = corr(first_maps{i}, second_maps{i}, 'Type', 'Kendall');
    Iself_real(i) = self_corr;

    others_corrs = zeros(N_participants - 1, 1);
    count = 1;
    for j = 1:N_participants
        if j ~= i
            others_corrs(count) = corr(first_maps{i}, second_maps{j}, 'Type', 'Kendall');
            count = count + 1;
        end
    end

    Identifiability_real(i) = self_corr - mean(others_corrs);
    Differentiability_real(i) = (self_corr - mean(others_corrs)) / std(others_corrs);
end


mean_Iself_real = mean(Iself_real);
mean_Identifiability_real = mean(Identifiability_real);
mean_Differentiability_real = mean(Differentiability_real);

fprintf('Baseline Identifiability: %.4f\n', mean_Identifiability_real);
fprintf('Baseline Differentiability: %.4f\n', mean_Differentiability_real);

% 2. Initialize impact per voxel 
impact_Identifiability = zeros(n_voxels,1);
impact_Differentiability = zeros(n_voxels,1);

% 2.5. Neighbors
[dim, xx, yy, zz, connmat, dtempl] = Omega_neighbors_ly(source);

% 3. Permutation per voxel
parfor v = 1:n_voxels
    % Maps
    second_maps_perm = second_maps;

    % Neighbors and central voxel
    neighbors = find(connmat(v,:) == 1);

    for n = 1:length(neighbors)
        voxel_idx = neighbors(n);
        voxel_values = zeros(N_participants, 1);
        for i = 1:N_participants
            voxel_values(i) = second_maps{i}(voxel_idx);
        end
        permuted_values = voxel_values(randperm(N_participants));
        for i = 1:N_participants
            second_maps_perm{i}(voxel_idx) = permuted_values(i);
        end
    end

    % Metrics post-permutation
    Iself_perm = zeros(N_participants,1);
    Identifiability_perm = zeros(N_participants,1);
    Differentiability_perm = zeros(N_participants,1);

    for i = 1:N_participants
        self_corr = corr(first_maps{i}, second_maps_perm{i}, 'Type', 'Kendall');
        Iself_perm(i) = self_corr;

        others_corrs = zeros(N_participants - 1, 1);
        count = 1;
        for j = 1:N_participants
            if j ~= i
                others_corrs(count) = corr(first_maps{i}, second_maps_perm{j}, 'Type', 'Kendall');
                count = count + 1;
            end
        end

        Identifiability_perm(i) = self_corr - mean(others_corrs);
        Differentiability_perm(i) = (self_corr - mean(others_corrs)) / std(others_corrs);
    end

    % Save differences with baseline
    impact_Identifiability(v) = mean(Identifiability_perm) - mean_Identifiability_real;
    impact_Differentiability(v) = mean(Differentiability_perm) - mean_Differentiability_real;

    if mod(v,100) == 0
        fprintf('Voxel %d/%d processed...\n', v, n_voxels);
    end
end

fprintf('Permutations finished .\n');

% 4. Visualization
figure;
subplot(2,1,1);
scatter(1:n_voxels, impact_Identifiability, 20, 'filled');
xlabel('Voxels');
ylabel('Δ Identifiability');
title('Voxel impact on Identifability');
grid on;

subplot(2,1,2);
scatter(1:n_voxels, impact_Differentiability, 20, 'filled');
xlabel('Voxels');
ylabel('Δ Differentiability (z-score)');
title('Voxel impact on Differentiability');
grid on;


%%
voxel_inside = find(source.inside==1);

source2 = source;
source2.avg.pow(voxel_inside) = impact_Identifiability;
source2.avg.mom=cell(length(source2.avg.noise),1);
source2.time=1;



cfg = [];
cfg.parameter  = 'avg.pow';
cfg.downsample = 2;
cfg.interpmethod  =  'nearest';
source_interp  = ft_sourceinterpolate (cfg, source2, source_forward.mri);
   

    figure('WindowState','maximized','Color',[1 1 1]);
    % figure

    cfg               = [];
    cfg.figure        = 'gca';
    cfg.method        = 'surface';
    cfg.funparameter  = 'pow';
    cfg.maskparameter = cfg.funparameter;
    cfg.funcolormap   = 'hot';
    % cfg.funcolorlim   = [-0.007 0];  %  scale idd


    cfg.projmethod    = 'nearest';
    cfg.opacity       = 0.8;
    cfg.camlight      = 'no';
    cfg.colorbar      = 'yes';
    cfg.surffile     = 'surface_pial_left.mat';
    cfg.surfinflated  = 'surface_inflated_left_caret_white.mat';
    subplot(2,2,1), ft_sourceplot(cfg,source_interp), view([-90 0]), camlight('left')
     

    subplot(2,2,3), ft_sourceplot(cfg,source_interp), view([90 0]),  camlight('left')



    cfg.surffile     = 'surface_pial_right.mat';
    cfg.surfinflated  = 'surface_inflated_right_caret_white.mat';
    subplot(2,2,2), ft_sourceplot(cfg,source_interp), view([90 0]),  camlight('right')

    subplot(2,2,4), ft_sourceplot(cfg,source_interp), view([-90 0]), camlight('right')



