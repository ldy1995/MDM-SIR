%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% demo of the iterative reconstraction method with E-3DTV prior in the
% following paper:
%   @ARTICLE{9115052,
%   author={Li, Danyang and Zeng, Dong and Li, Sui and Ge, Yongshuai and Bian, Zhaoying and Huang, Jing and Ma, Jianhua},
%   journal={IEEE Transactions on Medical Imaging},
%   title={MDM-PCCT: Multiple Dynamic Modulations for High-Performance Spectral PCCT Imaging},
%   year={2020},
%   volume={39},
%   number={11},
%   pages={3630-3642},
%   doi={10.1109/TMI.2020.3001616}}
%
% If you are interested in the method or have any queastions, please feel free
% to contact with Danayng Li (email: lidanyang1995@smu.edu.cn or dyli0730@gmail.com)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

%% load data
read_path = 'data/';
load([read_path, 'real_data_PCD.mat'])
addpath('toolbox/')

%% system matrix
SourceDetectorDistance = 1200;	    % in millimeters
SourceAxisDistance     = 850;	    % in millimeters
ProjectionNumber       = 480;
DetectorType           = inf;		% 0:Cylinder; inf:Flatpanel
DetectorTotalWidth     = 510.000;	% in millimeters
DetectorTotalHeight    = 6.000;		% in millimeters
DetectorColumnNumber   = 850;
DetectorRowNumber      = 10;
DetectorOffsetVert     = 0;			% in pixels
DetectorOffsetHoriz    = -111/3;    % in pixels
PixelNumberX           = 512; 	   % number of recon. slices in TOMO
PixelNumberY           = 512; 	   % In-plane pixel number

sg = sino_geom('fan', ...
    'nb', DetectorColumnNumber ,...
    'na', ProjectionNumber, ...
    'ds', 0.65,...%0.65
    'offset_s', DetectorOffsetHoriz, ... % channal offset in pixels
    'dsd', SourceDetectorDistance,...
    'dso', SourceAxisDistance,...
    'dfs', 0);  %arc 0   flat inf

ig = image_geom('nx', PixelNumberX,...
    'ny', PixelNumberY,...
    'fov', 150);

G = Gtomo2_dscmex(sg, ig);

%% parameter setting
% reconstruction parameters
par.alpha       = 0.3; % default [0.1, 1]
par.recon_iter  = 15; % default [10, 20]
par.v_min       = min(img_noise(:));
par.v_max       = max(img_noise(:));

% par for E-3DTV
par.rank        = [2,2,2];
par.rho         = 1.5;
par.maxIter     = 15;
par.mu          = 0.5;
par.lambda      = 0.4; % [0.3, 0.5] more smaller, more smoother
par.tau         = 0.5; % default [0.4, 0.6]

%% reconstrution
x0   = img_noise;
sino = sino_noise_PCD;

tic
img_E3DTV_denoise = spectral_E3DTV(sino, G, x0, par);
toc

% imshow
figure(100),
imshow([img_noise(:,:,1), img_noise(:,:,2), img_noise(:,:,3), ...
    img_noise(:,:,4);...
    img_E3DTV_denoise(:,:,1),img_E3DTV_denoise(:,:,2),img_E3DTV_denoise(:,:,3),...
    img_E3DTV_denoise(:,:,4)], [0.01, 0.03])
title(['recon-alpha-', num2str(par.alpha),...
    '-tau-', num2str(par.tau)])
fr = getframe(gcf);
I  = frame2im(fr);

% print the image
imwrite(I, ['img/Results_recon_alpha_', num2str(par.alpha),...
    '_tau_', num2str(par.tau),'.tif'], 'tif')
