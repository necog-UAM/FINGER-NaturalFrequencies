%% Beamforming of within-session group in Arana et al. (2025).

% This script is employed for source reconstruction of data from The Open MEG Archive. 
% Specifically for within-session group, the entire recording of each subject's continuous data is split in two halves before the source reconstruction.
% The output is both halves of the dataclean (dataclean1 and dataclean2) to treat them separately in the following steps and the repective sources forward and inverse for each half of the recording (source_forward_10mm_1 and source_forward_10mm_2).
% The input is mri_coreg and dataclean.

clear all
close all
clc

restoredefaultpath
addpath ('Z:\Toolbox\fieldtrip-20230118');
ft_defaults
addpath(genpath('Z:\Fingerprinting\scripts\Final'));


dpath = 'Z:\OMEGA\OMEGA_data\';
outpath = 'G:\Fingerprinting\Omega_data\'


% One sess, within-session group
subs = {'sub-0001' 'sub-0002' 'sub-0003' 'sub-0004' 'sub-0005' 'sub-0006' 'sub-0007' 'sub-0008' 'sub-0009' 'sub-0011' 'sub-0012' 'sub-0014' 'sub-0015' 'sub-0016' 'sub-0018' 'sub-0019' 'sub-0020' 'sub-0021' 'sub-0022' 'sub-0023' 'sub-0024' 'sub-0025' 'sub-0026' 'sub-0027' 'sub-0028' 'sub-0029' 'sub-0030' 'sub-0031' 'sub-0032' 'sub-0033' 'sub-0034' 'sub-0035' 'sub-0037' 'sub-0039' 'sub-0040' 'sub-0041' 'sub-0042' 'sub-0044' 'sub-0045' 'sub-0046' 'sub-0047' 'sub-0048' 'sub-0049' 'sub-0050' 'sub-0051' 'sub-0052' 'sub-0055' 'sub-0056' 'sub-0057' 'sub-0058' 'sub-0059' 'sub-0060' 'sub-0061' 'sub-0062' 'sub-0063' 'sub-0064' 'sub-0065' 'sub-0067' 'sub-0068' 'sub-0069' 'sub-0070' 'sub-0071' 'sub-0072' 'sub-0073' 'sub-0074' 'sub-0075' 'sub-0076' 'sub-0077' 'sub-0078' 'sub-0079' 'sub-0080' 'sub-0084' 'sub-0085' 'sub-0087' 'sub-0088' 'sub-0089' 'sub-0090' 'sub-0091' 'sub-0092' 'sub-0094' 'sub-0095' 'sub-0096' 'sub-0097' 'sub-0098' 'sub-0099' 'sub-0101' 'sub-0102' 'sub-0103' 'sub-0104' 'sub-0105' 'sub-0106' 'sub-0134' 'sub-0145' 'sub-0146' 'sub-0148' 'sub-0149' 'sub-0150' 'sub-0151' 'sub-0152' 'sub-0154' 'sub-0155' 'sub-0156' 'sub-0157' 'sub-0158' 'sub-0159' 'sub-0160' 'sub-0161' 'sub-0165' 'sub-0166' 'sub-0167' 'sub-0168' 'sub-0169' 'sub-0170' 'sub-0171' 'sub-0175' 'sub-0176' 'sub-0177' 'sub-0179' 'sub-0181' 'sub-0184' 'sub-0185' 'sub-0195' 'sub-0197' 'sub-0200' 'sub-0207' 'sub-0208' 'sub-0210' 'sub-0212'};
sess = {'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0002' 'ses-0003' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0003' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0002' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001' 'ses-0001'};


% 3.1. MRI normalization (output: mri_norm)
% 3.2. Head model
% 3.3. Dataclean split in two parts (output: dataclean1 and dataclean2) and forward model (output: source_forward_10mm_1 and source_forward_10mm_2)
% 3.4. Source reconstruction (output: source_inverse_10mm_1 and source_inverse_10mm_2)

for sub=1:length(subs)

    %% 3.1. MRI normalization
    % mri normalized (mrin) and transformation matrix (normtrans)

    cd([dpath  subs{sub} '\' sess{sub} ])
    load mri_coreg

    cfg            = [];
    cfg.nonlinear  = 'no';
    % cfg.spmversion = 'spm8';
    mrin           = ft_volumenormalise(cfg, mri);      % do you want to change the anatomical labels for the axes [Y, n]? Y (r,a,s,i)

    % determine the affine source->template coordinate transformation (from fieldtrip-20180405)
    normtrans = mrin.params.VG.mat * inv(mrin.params.Affine) * inv(mrin.params.VF.mat) * mrin.initial;

    mkdir([outpath  subs{sub} '\' sess{sub} ]);
    cd([outpath  subs{sub} '\' sess{sub} ])
    save mri_norm mrin normtrans

    %% 3.2. Head model
    % semi-realistic singleshell head model based on the implementation from Guido Nolte

    cfg             = [];
    cfg.spmversion  = 'spm8';
    segment         = ft_volumesegment(cfg,mri);      % extract brain surface
    segment.anatomy = mri.anatomy;

    %  Check that the segmentation is coregistered with mri
    figure
    cfg = [];
    cfg.interactive = 'yes';
    ft_sourceplot(cfg,mrin);      % only mri
    cfg.funparameter = 'gray';
    ft_sourceplot(cfg,segment);  % segmented gray matter on top

    pause;  % wait for user to close the figure


    cfg        = [];
    cfg.method = 'singleshell';
    vol        = ft_prepare_headmodel(cfg, segment);    % construct semi-realistic singleshell head model

    %% 3.3. Dataclean split in two parts and forward model
    % output: source_forward_10mm (contains mri, vol, grad, and the normalized grid and leadfields)

    cd([dpath  subs{sub} '\' sess{sub} ]);
    load dataclean

    half_length = round(length(dataclean.time{1,1}) / 2);       % divides time in half

    dataclean1 = dataclean;
    dataclean1.time = {dataclean.time{1,1}(1:half_length)};         % save first half 
    dataclean1.trial = {dataclean.trial{1,1}(:,1:half_length)};

    dataclean2 = dataclean;
    dataclean2.time = {dataclean.time{1,1}((half_length + 1):end)};       % save second half
    dataclean2.trial = {dataclean.trial{1,1}(:,(half_length + 1):end)};

    datacleans = {dataclean1, dataclean2};

    cd([outpath  subs{sub} '\' sess{sub} ]);
    save dataclean1
    save dataclean2

    cd([dpath  subs{sub} '\' sess{sub} ]);

    for dcleans =1:numel(datacleans)

        grad  = datacleans{dcleans}.grad;

        % Load normalized template grid (10mm)
        load standard_sourcemodel3d10mm
        grid = sourcemodel;

        % Load normalized mri and head model
        load mri_norm


        % Adapt the normalized grid to each individual's brain space
        posmni    = grid.pos;
        pos       = ft_warp_apply(inv(normtrans), grid.pos*10, 'homogenous')/10;
        grid.pos  = pos;
        grid.unit = 'cm';


        % % Convert grad, vol and grid to common units (mm)
        grad = ft_convert_units(grad, vol.unit);
        grid = ft_convert_units(grid, vol.unit);

        % Select only voxels within cortical mask (e.g. cerebellum is excluded) and corrected (inside the cortical surface projected with ft_sourceplot)
        % created with select_corticalvox_aal
        load ('correccion_vox_inside_10mm.mat')
        grid.inside = inside;

        % Compute leadfields for each grid's voxel
        cfg             = [];
        cfg.grid        = grid;
        cfg.grad        = grad;
        cfg.vol         = vol;
        cfg.channel     = {'MEG'};
        cfg.normalize   = 'yes';
        cfg.reducerank  = 2;
        grid2           = ft_prepare_leadfield(cfg);

        % Check that grad, vol and grid are correct (only for the first subject)

        figure
        plot3 (grad.chanpos(:,1), grad.chanpos(:,2), grad.chanpos(:,3), '.','MarkerEdgeColor',[0.8 0 0],'MarkerSize',25), hold on
        plot3 (vol.bnd.pos(:,1), vol.bnd.pos(:,2), vol.bnd.pos(:,3), '.','MarkerEdgeColor',[0 0 0.8]), hold on
        plot3 (grid2.pos(grid2.inside,1), grid2.pos(grid2.inside,2), grid2.pos(grid2.inside,3), '+k')

        pause;  % wait for user to close the figure

        % Save grad, vol, grid and mri in source_forward structure to be used later
        source_forward      = [];
        source_forward.vol  = vol;
        source_forward.mri  = mrin;
        source_forward.grad = grad;
        source_forward.grid = grid2;

        mkdir([outpath  subs{sub} '\' sess{sub} ]);
        cd([outpath  subs{sub} '\' sess{sub} ]);
        save(['source_forward_10mm_' num2str(dcleans) '.mat'], 'source_forward');


        %% 3.4. Computation of beamforming weights
        % output: source_inverse_10mm (contains beamforming weights in source.avg.filter)

        cfg            = [];
        cfg.covariance = 'yes';
        datacov        = ft_timelockanalysis(cfg, datacleans{dcleans});       % covariance matrix

        % Compute spatial filters (in source.avg.filter)
        cfg                   = [];
        cfg.method            = 'lcmv';
        cfg.grad              = source_forward.grad;
        cfg.headmodel         = source_forward.vol;
        cfg.grid              = source_forward.grid;
        cfg.lcmv.fixedori     = 'yes';
        cfg.lcmv.normalize    = 'yes';
        cfg.lcmv.projectnoise = 'yes';
        cfg.lcmv.keepfilter   = 'yes';          % important: save filters to use them later
        cfg.lcmv.lambda       = '10%';          % the higher the smoother
        cfg.lcmv.reducerank   = 2;
        source                = ft_sourceanalysis(cfg, datacov);

        load standard_sourcemodel3d10mm
        source.avg.ori = {};
        source.avg.mom = {};
        source.avg.noisecov = {};
        source.pos     = sourcemodel.pos;            % standard grid positions
        source.inside  = grid.inside;

        cd([outpath  subs{sub} '\' sess{sub} ]);
        save(['source_inverse_10mm_' num2str(dcleans) '.mat'], 'source');        
    end
end
