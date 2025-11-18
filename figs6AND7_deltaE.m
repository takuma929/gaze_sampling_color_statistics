%% ========================================================================
%  Script: figs6AND7_deltaE.m
%
%  Purpose:
%    Generate ΔE plots for human observers and computational observers
%    across:
%        - 7 selected scenes
%        - 50-scene extended set
%
%    Includes:
%       • Loading of ΔE (CIECAM16-UCS) adaptation summary
%       • Extraction of human & model data
%       • Ordering computational observers by performance
%       • Creating Fig6 (7 scenes) and Fig7 (50 scenes)
%       • Generating montage thumbnails of natural scenes
%
%  Inputs:
%       data/DeltaE_adaptation_summary.mat   → de_ciecam16UCS_median
%       data/allData.mat                        → raw images
%       data/fig_parameters.mat              → typography, layout
%       data/list_Obs_order.mat              → mapping of simulation IDs
%
%  Outputs:
%       /figs/Fig6_tau1.pdf
%       /figs/Fig6_tau10.pdf
%       /figs/Fig7_tau1.pdf
%       /figs/Fig7_tau10.pdf
%       /figs/Fig7_thumbnail.png
%
%  ========================================================================

clearvars; close all; clc;


%% ========================================================================
% Load Required Data
% ========================================================================
load(fullfile(pwd,'data','DeltaE_adaptation_summary.mat'), ...
     'de_ciecam16UCS_median');

load(fullfile(pwd,'data','allData.mat'), 'allData');
load(fullfile(pwd,'data','fig_parameters.mat'), 'figp');
load(fullfile(pwd,'data','list_Obs_order'), 'list_Obs');

%% ========================================================================
% Scene Lists
% ========================================================================
sceneList_all = fieldnames(de_ciecam16UCS_median);

% First 7 scenes are "7-scenes set"
sceneList = sceneList_all(1:7);

% Excluded indices from the 50-scene list
sceneN_exclude = [8 9 14 40 50];

% Build the clean 50-scene list (skipping excluded)
sceneList_50scenes_temp = sceneList_all(8:end);
sceneList_50scenes = {};
cnt = 0;
for N = 1:50
    if ~ismember(N, sceneN_exclude)
        cnt = cnt + 1;
        sceneList_50scenes{cnt} = sceneList_50scenes_temp{N};
    end
end


%% ========================================================================
% Observer Model Definition
% ========================================================================
list_compObs = { ...
    'randomWalk_step_onethirdmode', ...
    'randomWalk_step_mode', ...
    'randomWalk_step_3xmode', ...
    'randomGaze', ...
    'maxChroma'};

% Just-noticeable difference threshold (CIECAM16-UCS)
jnd = 2.46;


%% ========================================================================
% Extract ΔE for 7 scenes (humans + computational observers)
% ========================================================================
data_human_mean = zeros(6,2,7);
data_compObs = struct();

sceneN = 0;
for s = sceneList'
    sceneN = sceneN + 1;
    name = s{1};

    % -- Average human across observers
    data_human_mean(:,:,sceneN) = mean(de_ciecam16UCS_median.(name).cct4000.human,3);

    % -- Simulations
    data_sim = de_ciecam16UCS_median.(name).cct4000.simulation;

    for compObs = list_compObs
        model = compObs{1};
        id = find(contains(list_Obs.simulation, model));
        data_compObs.(model)(:,:,sceneN) = mean(data_sim(:,:,id),3);
    end
end

% Std across scenes for humans
data_human_std = squeeze(std(data_human_mean,[],3));


%% ========================================================================
% Extract ΔE for 50 scenes (computational observers only)
% ========================================================================
sceneN = 0;
data_compObs_50scenes = struct();

for s = sceneList_50scenes
    sceneN = sceneN + 1;
    name = s{1};

    data_sim = de_ciecam16UCS_median.(name).cct4000.simulation;

    for compObs = list_compObs
        model = compObs{1};
        id = find(contains(list_Obs.simulation, model));
        data_compObs_50scenes.(model)(:,:,sceneN) = mean(data_sim(:,:,id),3);
    end
end


%% ========================================================================
% Frame/Time Axes (fixed, preserved from original script)
% ========================================================================
frameN_max = 14999;
ndata = 6;
list_frameN = round(10.^(linspace(0,log10(frameN_max),ndata)));
list_frameN(1) = 3;   % use 3rd frame instead of first
list_time = list_frameN / 50;


%% ========================================================================
% Plot Styling (unchanged)
% ========================================================================
cmap.randomWalk_step_onethirdmode = [0.7412 0.8431 0.9059];
cmap.randomWalk_step_mode         = [0.4196 0.6824 0.8392];
cmap.randomWalk_step_3xmode       = [0.1294 0.4431 0.7098];
cmap.randomGaze                   = 'r';
cmap.maxChroma                    = [25 107 36]/255;
cmap.human                        = 'm';

symbol.randomGaze                   = 's'; markersize.randomGaze = 6;
symbol.randomWalk_step_onethirdmode = 'd'; markersize.randomWalk_step_onethirdmode = 6;
symbol.randomWalk_step_mode         = '^'; markersize.randomWalk_step_mode = 6;
symbol.randomWalk_step_3xmode       = '<'; markersize.randomWalk_step_3xmode = 6;
symbol.maxChroma                    = '>'; markersize.maxChroma = 6;

tau_label = {'tau1','tau10'};


%% ========================================================================
% Main Loop: Plot for 7 scenes and 50 scenes
% ========================================================================
for sceneName = {'7scenes','50scenes'}

    disp(['Plotting for... ' sceneName{1}]);

    for tau = 1:2

        fig = figure; hold on;

        % Identify scene index
        if strcmp(sceneName{1},'7scenes')
            idxScene = strcmp(sceneList, sceneName{1});
        else
            idxScene = 8; % index for 50-scenes (placeholder in original)
        end

        % ------------------------------------------------------
        % Plot computational observers
        % ------------------------------------------------------
        for compObs = list_compObs
            model = compObs{1};

            if strcmp(sceneName{1},'50scenes')
                y = mean(data_compObs_50scenes.(model)(:,tau,:),3);
                se = std(data_compObs_50scenes.(model)(:,tau,:),[],3)/sqrt(45);
            else
                y = mean(data_compObs.(model)(:,tau,:),3);
                se = [];
            end

            x = log10(list_time);
            plot(x, y, ...
                 strcat(symbol.(model),'-'), ...
                 'MarkerSize', markersize.(model), ...
                 'MarkerEdgeColor', ones(3,1)*0.97, ...
                 'MarkerFaceColor', cmap.(model), ...
                 'Color', cmap.(model), ...
                 'LineWidth', 0.5);

            if strcmp(sceneName{1},'50scenes')
                e = errorbar(x,y,se,se); 
                e.CapSize = 0; e.LineStyle = 'none'; e.Color = cmap.(model);

                temp = squeeze(data_compObs_50scenes.(model)(:,tau,:));
                [sortedDeltaE.(model), rank.(model)] = sort(min(temp));
            end
        end

        % ------------------------------------------------------
        % Plot human (7 scenes only)
        % ------------------------------------------------------
        if strcmp(sceneName{1},'7scenes')
            x = log10(list_time);
            Y = data_human_mean(:,tau,:);
            y = mean(Y,3);
            se = std(Y,[],3)/sqrt(7);

            plot(x,y,'o-','LineWidth',1,'MarkerSize',7, ...
                 'MarkerEdgeColor',ones(3,1)*0.97,'MarkerFaceColor','k','Color','k');

            e = errorbar(x,y,se,se);
            e.CapSize = 0; e.LineStyle = 'none'; e.Color = 'k';

            minDeltaE_5deg = squeeze(min(data_human_mean))';
        end

        %% ------------------------------------------------------
        %% Axes + Figure Formatting  (EXACTLY preserved)
        %% ------------------------------------------------------
        ax = gca;
        ax.XTick = [-2 -1 0 1 2]; ax.XLim = [-1.5 2.5];
        ax.XTickLabel = char('-2','-1','0','1','2');

        ax.YTick = [0 5.0 10]; ax.YLim = [0 11];
        ax.YTickLabel = char('0.0','5.0','10.0');

        fig.PaperType      = 'a4';
        fig.PaperUnits     = 'centimeters';
        fig.Units          = 'centimeters';
        fig.Color          = 'w';
        fig.InvertHardcopy = 'off';
        fig.PaperPosition  = [0,10,8.45,8.45];
        fig.Position       = [10,10,figp.twocolumn/4*0.95,figp.twocolumn/4*0.95];

        ax.FontName = figp.fontname;
        ax.FontSize = figp.fontsize;
        ax.Color    = ones(3,1)*0.97;

        xlabel('', 'FontSize', figp.fontsize_axis);
        ylabel('', 'FontSize', figp.fontsize_axis);

        % jnd reference lines
        for level = 1:4
            r = refline(0, jnd*level);
            r.Color = 'm'; r.LineStyle = ':'; r.LineWidth = 0.5;
        end

        ax.Units = 'centimeters';
        axis square;
        ax.Position = [0.65 0.65 3.2 3.2];
        ax.TickDir = 'out';
        ticklengthcm(ax, 0.1);
        box off; grid off;
        ax.LineWidth = 0.3;

        %% ------------------------------------------------------
        %% Save Figures
        %% ------------------------------------------------------
        if strcmp(sceneName{1},'7scenes')
            fname = ['Fig6_',tau_label{tau}];
        else
            fname = ['Fig7_',tau_label{tau}];
        end

        exportgraphics(fig, fullfile('figs',[fname,'.pdf']), ...
            'ContentType','vector');
    end

    close all
end


%% ========================================================================
% Make Thumbnail for Fig 7
% ========================================================================
disp('Making thumbnail...');

imgs_all = [];
imgs_horz = [];
cnt = 0;

% Order scenes by rank of randomWalk_step_3xmode
orderedScenes = sceneList_50scenes(rank.randomWalk_step_3xmode);

for sceneName = orderedScenes
    cnt = cnt + 1;
    img = allData.(sceneName{1}).bmp;

    img(:,:,1) = img(:,:,1)*0.9;    % whiten balance adjustment
    imgs_horz = [imgs_horz img];

    if mod(cnt,10)==0
        imgs_all = [imgs_all; imgs_horz];
        imgs_horz = [];
    end

    if cnt==45 && ~(mod(cnt,10)==0)
        padSize = size(img,2)*(10-mod(cnt,10));
        imgs_all = [imgs_all; ...
            [imgs_horz, 256*ones(size(imgs_horz,1), padSize, 3)]];
    end
end

thumbnail = imresize(imgs_all * 1.1, 0.25);
imwrite(thumbnail, fullfile('figs','Fig7_thumbnail.png'));

disp('All figures saved successfully.');
