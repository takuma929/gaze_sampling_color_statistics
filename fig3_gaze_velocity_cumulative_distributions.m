%% ========================================================================
%  Script: fig3_gaze_velocity_cumulative_distributions.m
%
%  Purpose:
%    For each observer (FIELD and LAB measurements), compute and plot the
%    cumulative distribution of gaze velocities, averaged across all 7
%    scenes. The script:
%
%      (1) Loads gaze trajectories for each observer in each scene
%      (2) Computes instantaneous gaze velocity (deg/s)
%      (3) Builds a histogram and cumulative distribution
%      (4) Averages across scenes, computes SEM (or SD here)
%      (5) Plots mean ± SD shading for each observer
%      (6) Saves the plots as vector PDF (publication-ready)
%
%  Inputs:
%    allData (MAT-file) – contains scene images + gaze data
%    figp (MAT-file)    – figure style parameters (font, size, width)
%
%  Outputs:
%    ./fig/Fig3_speed_histogram_field.pdf
%    ./fig/Fig3_speed_histogram_lab.pdf
%
%  Notes:
%    • Velocities are scaled by factor 50 / 48.6957 (deg/unit → deg/s)
%    • Cumulative percentages computed per scene, then averaged
%
%  Author: Takuma Morimoto
%  Updated: 18th Nov 2025
% ========================================================================

clearvars; close all; clc;


%% ------------------------------------------------------------------------
% Set repository path and load data
% -------------------------------------------------------------------------
load(fullfile(pwd,'data','allData'),'allData');   % main data struct
load(fullfile(pwd,'data','fig_parameters'),'figp');   % figure parameters


%% ------------------------------------------------------------------------
% Scene & subject definitions
% -------------------------------------------------------------------------
sceneList = fieldnames(allData);       % should contain 7 scenes

subjectList.field = {'dhf','ka'};
subjectList.lab   = {'dhf','ka','eo','ma','nm'};

%% ------------------------------------------------------------------------
% Colormap
% -------------------------------------------------------------------------
cmap = [0.3600,0.6847,0.5824;
    0.8153,0.4871,0.6882;
    0.8894,0.4976,0.3459;
    0.4976,0.5647,0.7165;
    0.5859,0.7624,0.2965];

%% ------------------------------------------------------------------------
% Speed histogram parameters
% -------------------------------------------------------------------------
numBins  = 1000;
maxSpeed = 20;                          % unify speed range
binEdges = linspace(0, maxSpeed, numBins + 1);
binCenters = binEdges(1:end-1) + diff(binEdges)/2;

%% ========================================================================
% Main loop: FIELD and LAB
% ========================================================================
for env = {'field','lab'}
    currentEnv = env{1};

    % Create figure
    fig = figure; hold on;

    subjectCounter = 0;
    maxSpeedCounter = 1;

    % Loop through subjects
    for subjectCell = subjectList.(currentEnv)
        subjectName = subjectCell{1};
        subjectCounter = subjectCounter + 1;

        % Storage for cumulative distributions across 7 scenes
        cumCountsAcrossScenes = zeros(7, numBins);

        % -----------------------------------------------------------------
        % Loop through scenes
        % -----------------------------------------------------------------
        for s = 1:7
            sceneName = sceneList{s};

            % Extract gaze positions
            x = allData.(sceneName).gazeData.human.(currentEnv).(subjectName).xdata;
            y = allData.(sceneName).gazeData.human.(currentEnv).(subjectName).ydata;

            % Compute velocities
            dx = diff(x);
            dy = diff(y);
            speed = sqrt(dx.^2 + dy.^2) * (50 / 48.6957);   % convert to deg/s

            % Store maximum speed
            maxSpeedObserved.(currentEnv)(maxSpeedCounter,1) = max(speed);
            maxSpeedCounter = maxSpeedCounter + 1;

            % Histogram + cumulative (%)
            counts = histcounts(speed, binEdges);
            cumCountsAcrossScenes(s,:) = cumsum(counts) / length(speed) * 100;
        end

        % -----------------------------------------------------------------
        % Compute mean ± SD across scenes
        % -----------------------------------------------------------------
        meanCum = mean(cumCountsAcrossScenes, 1);
        sdCum   = std(cumCountsAcrossScenes, 0, 1);   % SD across scenes

        % Prepare shading polygon
        shadeX = [binCenters, fliplr(binCenters)];
        shadeY = [meanCum + sdCum, fliplr(meanCum - sdCum)];

        % -----------------------------------------------------------------
        % Plot shading + mean curve
        % -----------------------------------------------------------------
        fill(shadeX, shadeY, cmap(subjectCounter,:), ...
            'FaceAlpha', 0.30, 'EdgeColor', 'none');
        plot(binCenters, meanCum, ...
            'Color', cmap(subjectCounter,:), 'LineWidth', 1.2);

        % Print value at 5 deg/s (debug/info)
        [~, idx5deg] = min(abs(binCenters - 5));
        fprintf('%s (%s): %.1f%% at 5 deg/s\n', subjectName, currentEnv, meanCum(idx5deg));
    end


    %% --------------------------------------------------------------------
    % Axes formatting
    % ---------------------------------------------------------------------
    xlabel('Gaze velocity [deg s^{-1}]', 'FontSize', figp.fontsize_axis);
    ylabel('Cumulative frequency [%]',   'FontSize', figp.fontsize_axis);

    xlim([-0.5, maxSpeed]);
    ylim([0, 100]);

    % Reference vertical line at 5 deg/s
    line([5, 5], [0 150], 'Color', 'm', ...
         'LineStyle', ':', 'LineWidth', 0.5);

    ax = gca;
    ax.XTick = 0:10:maxSpeed;
    ax.YTick = 0:25:100;
    ax.FontName = figp.fontname;
    ax.FontSize = figp.fontsize;
    ax.TickDir = 'in';
    ax.LineWidth = 0.25;
    box off;
    grid on; ax.XGrid = 'off';
    ax.XColor = 'k'; ax.YColor = 'k';


    %% --------------------------------------------------------------------
    % Figure sizing and export
    % ---------------------------------------------------------------------
    fig.PaperType  = 'a4';
    fig.PaperUnits = 'centimeters';
    fig.Units      = 'centimeters';
    fig.Color      = 'w';
    fig.Position   = [10, 10, figp.twocolumn/2, figp.twocolumn/3];

    outFile = sprintf('Fig3_speed_histogram_%s.pdf', currentEnv);
    exportgraphics(fig, fullfile('figs', outFile), ...
        'ContentType', 'vector', 'BackgroundColor', 'none');
end

disp('All histograms (Fig 3) generated and saved.');
