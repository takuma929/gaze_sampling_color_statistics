%% ========================================================================
%  Script: figs1AND2_gaze_density.m
%
%  Purpose:
%    Generate gaze-density montage figures for each natural scene in both
%    the FIELD and LAB measurements. For each observer, the script:
%       (1) loads the scene image
%       (2) gamma-corrects and rescales it for optimal visibility
%       (3) overlays gaze positions as semi-transparent white points
%       (4) assembles a horizontal montage of the base image + gaze plot
%    The final montages are saved as publication-ready PNG files.
%
%  Inputs (loaded automatically):
%    allData  – structure containing scene images and gaze coordinates
%    figp     – structure containing figure parameters (unused but loaded)
%
%  Output:
%    PNG files saved to ./figs/
%      Fig1_gazeDensity_[scene].png   (FIELD environment)
%      Fig2_gazeDensity_[scene].png   (LAB environment)
%
%  Author: Takuma Morimoto
%  Last updated: 18th Nov 2025
% =========================================================================

clearvars; close all; clc;     % Clean workspace and command window

%% ------------------------------------------------------------------------
% Load required data files
% -------------------------------------------------------------------------
load(fullfile(pwd,'data','allData.mat'),'allData');   % scene images + gaze data
load(fullfile(pwd,'data','fig_parameters.mat'),'figp'); % general fig params

%% ------------------------------------------------------------------------
% Scene list and environment/observer definitions
% -------------------------------------------------------------------------
sceneList = {'cp3','elevator','forest','gap','orange','ruivaes','trunks'};
envList   = {'field','lab'};

observerList.field = {'dhf','ka'};
observerList.lab   = {'dhf','ka','eo','ma','nm'};

%% ------------------------------------------------------------------------
% Gamma and scaling parameters for each scene
% (tuned manually for publication-quality contrast)
% -------------------------------------------------------------------------
scale.cp3     = 1.1;  gamma.cp3     = 1.0;
scale.elevator = 1.0; gamma.elevator = 1.0;
scale.forest   = 1.2; gamma.forest   = 0.6;
scale.gap      = 1.0; gamma.gap      = 1.0;
scale.orange   = 1.0; gamma.orange   = 1.0;
scale.ruivaes  = 1.0; gamma.ruivaes  = 0.8;
scale.trunks   = 1.0; gamma.trunks   = 1.0;

%% ------------------------------------------------------------------------
% Main loop: iterate over environments (FIELD / LAB) and scenes
% -------------------------------------------------------------------------
for env = envList
    currentEnv = env{1};

    for sceneCell = sceneList
        sceneName = sceneCell{1};

        fprintf('Processing scene: %s (%s environment)\n', sceneName, currentEnv);

        %% ---------------------------------------------------------------
        % Load and prepare base scene image
        % ---------------------------------------------------------------
        baseImage = double(allData.(sceneName).bmp) / 255;             % normalize
        baseImage = baseImage .^ gamma.(sceneName);                    % gamma correction
        baseImage = baseImage ./ max(baseImage(:));                    % renormalize
        baseImage = baseImage * scale.(sceneName);                     % scale for visibility

        % Initialize montage with the base image
        montageImage = baseImage;

        %% ---------------------------------------------------------------
        % Loop through all observers in this environment
        % ---------------------------------------------------------------
        for obs = observerList.(currentEnv)
            observerName = obs{1};

            % Extract gaze coordinates
            xGaze = allData.(sceneName).gazeData.human.(currentEnv).(observerName).xdata;
            yGaze = allData.(sceneName).gazeData.human.(currentEnv).(observerName).ydata;

            % Create figure for gaze overlay
            fig = figure('Visible','off');
            imshow(baseImage * 0.5); hold on;    % slightly dim background
            scatter(xGaze, yGaze, 10, 'wo', 'filled', ...
                'MarkerFaceAlpha', 0.10);        % semi-transparent gaze dots
            hold off;

            % Capture plot as an image
            gazeFrame = getframe(gca);
            gazeImage = double(gazeFrame.cdata) / 255;

            % Resize gaze image to match original scene height
            gazeImage = imresize(gazeImage, size(baseImage,1:2));

            % Add a white separator + gaze image to montage
            separator = ones(size(baseImage,1), 10, 3);
            montageImage = [montageImage separator gazeImage];

            close(fig);
        end

        %% ---------------------------------------------------------------
        % Save output figure
        % ---------------------------------------------------------------
        outDir = fullfile(pwd,'figs');
        if ~exist(outDir, 'dir')
            mkdir(outDir);
        end

        if strcmp(currentEnv,'lab')
            outName = sprintf('Fig2_gazeDensity_%s.png', sceneName);
        else
            outName = sprintf('Fig1_gazeDensity_%s.png', sceneName);
        end

        imwrite(montageImage, fullfile(outDir, outName));
    end
end

disp('All gaze-density figures (Figs 1 and 2) generated and saved.');
