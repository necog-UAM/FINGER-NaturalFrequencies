%% Beamforming of between-session group in Arana et al. (2025).

% This script is employed for source reconstruction of data from The Open MEG Archive. 
% The outputs are the source forward and inverse.
% The input is mri_coreg and dataclean.
% You should select the first session or the second.

clear all
close all
clc

restoredefaultpath
addpath ('Z:\Toolbox\fieldtrip-20230118');
ft_defaults
addpath(genpath('Z:\Fingerprinting\scripts\Final'));


dpath = 'Z:\OMEGA\OMEGA_data\';
rawpath = 'Z:\OMEGA\OMEGA_raw\';

% Subjects of between-session group that have 2 sessions
subs = {'sub-0001'	'sub-0002'	'sub-0006'	'sub-0008'	'sub-0011'	'sub-0016'	'sub-0019'	'sub-0020'	'sub-0022'	'sub-0023'	'sub-0025'	'sub-0030'	'sub-0032'	'sub-0035'	'sub-0039'	'sub-0040'	'sub-0041'	'sub-0042'	'sub-0044'	'sub-0046'	'sub-0048'	'sub-0049'	'sub-0050'	'sub-0051'	'sub-0106'	'sub-0150'	'sub-0200'};

% Select session
select_session = 1;   % 1 for fisrt session, 2 for second session

if select_session == 1

    % First sess
    sess = {'ses-0001' 'ses-0004' 'ses-0001'  'ses-0001'	'ses-0002'	'ses-0002'	'ses-0003'	'ses-0004'	'ses-0002'	'ses-0001'	'ses-0003'	'ses-0001'	'ses-0001'	'ses-0003'	'ses-0001'	'ses-0005'	'ses-0003'	'ses-0001'  'ses-0002'	'ses-0002'	'ses-0004'	'ses-0002'	'ses-0004'	'ses-0001'  'ses-0001'	'ses-0001'	'ses-0001'};

elseif select_session == 2

    % Second sess
    sess = {'ses-0003' 'ses-0002' 'ses-0003'	'ses-0002'	'ses-0001' 	'ses-0001'	'ses-0002'	'ses-0003'	'ses-0001'	'ses-0002'	'ses-0002'	'ses-0002'	'ses-0002'	'ses-0002'	'ses-0002'	'ses-0004'	'ses-0001'	'ses-0002'	'ses-0003'	'ses-0001'	'ses-0002'	'ses-0001'	'ses-0002'	'ses-0002'	'ses-0002'	'ses-0002'	'ses-0002'};

else
    disp('You should select session 1 or 2')
end


% 3.1. MRI normalization (output: mri_norm)
% 3.2. Head model
% 3.3. Forward model (output: source_forward_10mm)
% 3.4. Source reconstruction (output: source_inverse_10mm)

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


save mri_norm mrin normtrans

%% 3.2. Head model
% semi-realistic singleshell head model based on the implementation from Guido Nolte

cfg             = [];
cfg.spmversion  = 'spm8';
segment         = ft_volumesegment(cfg,mri);      % extract brain surface
segment.anatomy = mri.anatomy; 

% Check that the segmentation is coregistered with mri
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

%% 3.3. Forward model
% output: source_forward_10mm (contains mri, vol, grad, and the normalized grid and leadfields)  

cd([dpath  subs{sub} '\' sess{sub} ])
load dataclean
grad  = dataclean.grad;

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

% % Check that grad, vol and grid are correct (only for the first subject)

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

save source_forward_10mm source_forward

%% 3.4. Computation of beamforming weights
% output: source_inverse_10mm (contains beamforming weights in source.avg.filter) 

cfg            = [];
cfg.covariance = 'yes';
datacov        = ft_timelockanalysis(cfg, dataclean);       % covariance matrix

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

cd([dpath  subs{sub} '\' sess{sub} ])
save source_inverse_10mm source
end
