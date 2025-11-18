%% ========================================================================
%  Script: figs4AND5_corr_localANDglobal.m
%
%  Purpose:
%    Generate field (outdoor) and lab (indoor) plots of the correlation
%    between LOCAL vs. GLOBAL signals as a function of
%    time. For each subject and each cone
%    class (L/M/S), the script:
%
%       (1) Loads pre-computed local/global correlation values
%       (2) Averages across 7 scenes
%       (3) Computes standard errors
%       (4) Overlays computational observer predictions
%       (5) Exports publication-ready vector PDF figures
%
%  Inputs:
%    corrcoeff_localvsglobal  – struct of correlation values (loaded)
%    list_frameN              – vector of frame window sizes in samples
%    figp                     – structure of typography + layout settings
%
%  Outputs:
%    /figs/Fig5_logL.pdf, Fig5_logM.pdf, Fig5_logS.pdf (FIELD)
%    /figs/Fig4_logL.pdf, Fig4_logM.pdf, Fig4_logS.pdf (LAB)
%
%  ------------------------------------------------------------------------
%  Author: Takuma Morimoto
%  Last updated: 18th Nov 2025
%  ========================================================================

clearvars; close all; clc;

%% Load data
load(fullfile(pwd,'data','corrcoeff_localvsglobal_thre_5deg'), ...
     'corrcoeff_localvsglobal','list_frameN');

load(fullfile(pwd,'data','fig_parameters.mat'),'figp');


%% ------------------------------------------------------------------------
% Subject lists / variable lists
% ------------------------------------------------------------------------
subjectList.field = {'dhf','ka'};
subjectList.lab   = {'dhf','ka','eo','ma','nm'};

varList.field = {'logL','logM','logS'};
varList.lab   = {'logL','logM','logS'};

%% ------------------------------------------------------------------------
% Colors
% ------------------------------------------------------------------------
cmap.field = [
    0.3600 0.6847 0.5824
    0.8153 0.4871 0.6882
    0.8894 0.4976 0.3459
    0.4976 0.5647 0.7165
    0.5859 0.7624 0.2965];

cmap.lab = cmap.field;  % same palette

compObsList = {'randomGaze','randomWalk_step_mode','maxChroma'};
cmap_compObs.randomGaze          = [1 0 0];
cmap_compObs.randomWalk_step_mode = [0 0 1];
cmap_compObs.maxChroma           = [0.0980,0.4196,0.1412];

%% ------------------------------------------------------------------------
% Helper: load 5-iteration computational observer data
% ------------------------------------------------------------------------
function M = loadCompObsMean(corrStruct, env, varName, compName)
    for i = 1:5
        fieldName = sprintf('%s_iter%d', compName, i);
        tmp(:,:,i) = squeeze(corrStruct.(env).(fieldName).(varName)(:,3,:)); %#ok<AGROW>
    end
    M = mean(tmp,3);   % average across 5 iterations
end


%% ========================================================================
% ----------------------------- FIELD -------------------------------------
% ========================================================================
disp('Plotting FIELD...');

env = 'field';
frameX = log10(list_frameN/50)';

for v = varList.(env)
    varName = v{1};

    % Precompute computational observers
    for c = compObsList
        comp = c{1};
        compObsMean.(comp).(varName) = ...
            loadCompObsMean(corrcoeff_localvsglobal, env, varName, comp);
    end

    fig = figure;
    hold on;

    % Subjects
    for s = 1:numel(subjectList.(env))
        subj = subjectList.(env){s};

        % 7-scene responses
        M = squeeze(corrcoeff_localvsglobal.(env).(subj).(varName)(:,3,:));
        y  = mean(M,2);
        se = std(M,[],2) / sqrt(size(M,2));

        % Subject trace + patch (COLOR & LINEWIDTH unchanged)
        plot(frameX, y, 'Color', cmap.(env)(s,:), 'LineWidth', 0.5); hold on
        patch([frameX; flip(frameX)], [y-se; flip(y+se)], cmap.(env)(s,:), ...
              'FaceAlpha',0.1, 'EdgeColor','none');

        % Only the first subject overlays computational observers
        if s == 1
            for c = compObsList
                comp = c{1};
                plot(frameX, mean(compObsMean.(comp).(varName),2), ...
                     'Color', cmap_compObs.(comp), ...
                     'LineWidth', 1, 'LineStyle', ':');
            end
        end

        % --- White backed text (unchanged) ------------------------------
        drawTextOpaque(gca, 2.1, -s*0.1+0.24, ...
                       sprintf('%.2f', max(y)), cmap.(env)(s,:), 7, 'Arial');
    end

    % --- Axes formatting (EXACT original) --------------------------------
    ax = gca;
    ax.XTick = [0 1 2]; ax.XLim = [-0.1 2.45];
    ax.YTick = 0:0.25:1; ax.YLim = [-0.1 1.05];
    ax.XTickLabel = char('0.0','1.0','2.0','2.5');
    ax.YTickLabel = char('0.0','','0.5','','1.0');

    fig.PaperType = 'a4'; fig.PaperUnits = 'centimeters';
    fig.Units = 'centimeters'; fig.Color = 'w';
    fig.InvertHardcopy = 'off';
    fig.PaperPosition = [0,10,8.45,8.45];
    fig.Position = [10,10,figp.twocolumn/3*0.95,figp.twocolumn/3*0.95];

    ax.FontName = figp.fontname; ax.FontSize = figp.fontsize;
    ax.Color = ones(3,1)*1;
    xlabel('','FontSize',figp.fontsize_axis); ylabel('','FontSize',figp.fontsize_axis);
    ax.Units = 'centimeters';
    axis square; ax.Position = [0.82 0.82 4.2 4.2];
    ax.TickDir = 'in'; ticklengthcm(ax,0.1)
    box off; grid on; ax.XGrid = 'off';
    ax.LineWidth = 0.25; ax.XColor = 'k'; ax.YColor = 'k';

    exportgraphics(fig, fullfile('figs',['Fig5_',varName,'.pdf']), ...
                   'ContentType','vector','BackgroundColor','none')
end


%% ========================================================================
% ----------------------------- LAB ---------------------------------------
% ========================================================================
disp('Plotting LAB...');

env = 'lab';
frameX = log10(list_frameN/50)';

for v = varList.(env)
    varName = v{1};

    fig = figure;
    hold on;

    for s = 1:numel(subjectList.(env))
        subj = subjectList.(env){s};

        M = squeeze(corrcoeff_localvsglobal.(env).(subj).(varName)(:,3,:));
        y  = mean(M,2);
        se = std(M,[],2) / sqrt(size(M,2));

        plot(frameX, y, 'Color', cmap.(env)(s,:), 'LineWidth', 0.5); hold on
        patch([frameX;flip(frameX)], [y-se;flip(y+se)], cmap.(env)(s,:), ...
              'FaceAlpha',0.2, 'EdgeColor','none');

        if s == 1
            for c = compObsList
                comp = c{1};
                plot(frameX, mean(compObsMean.(comp).(varName),2), ...
                     'Color', cmap_compObs.(comp), ...
                     'LineWidth',1, 'LineStyle',':');
            end
        end

        drawTextOpaque(gca, 2.1, -s*0.1+0.54, ...
                       sprintf('%.2f',max(y)), cmap.(env)(s,:), 7,'Arial');
    end

    % --- Axes formatting (unchanged) ------------------------------------
    ax = gca;
    ax.XTick = [0 1 2]; ax.XLim = [-0.1 2.45];
    ax.YTick = 0:0.25:1; ax.YLim = [-0.1 1.05];
    ax.XTickLabel = char('0.0','1.0','2.0','2.5');
    ax.YTickLabel = char('0.0','','0.5','','1.0');

    fig.PaperType = 'a4'; fig.PaperUnits = 'centimeters';
    fig.Units = 'centimeters'; fig.Color = 'w';
    fig.InvertHardcopy = 'off';
    fig.PaperPosition = [0,10,8.45,8.45];
    fig.Position = [10,10,figp.twocolumn/3*0.95,figp.twocolumn/3*0.95];

    ax.FontName = figp.fontname; ax.FontSize = figp.fontsize;
    ax.Color = ones(3,1)*1;
    xlabel('','FontSize',figp.fontsize_axis); ylabel('','FontSize',figp.fontsize_axis);
    ax.Units = 'centimeters'; axis square;
    ax.Position = [0.82 0.82 4.2 4.2];
    ax.TickDir = 'in'; ticklengthcm(ax,0.1)
    box off; grid on; ax.XGrid = 'off';
    ax.LineWidth = 0.25; ax.XColor = 'k'; ax.YColor = 'k';

    exportgraphics(fig, fullfile('figs',['Fig4_',varName,'.pdf']), ...
                   'ContentType','vector','BackgroundColor','none')
end

close all


%% ========================================================================
% Helper function
% ========================================================================
function drawTextOpaque(ax, x, y, str, color, fs, fontname)
    t = text(ax, x, y, str, 'FontSize', fs, 'FontName', fontname, ...
             'Color', color, 'Units','data', 'Visible','off');
    drawnow limitrate
    e = get(t,'Extent'); delete(t)
    padX = 0.04; padY = 0.025;
    rectangle(ax, 'Position', [e(1)-padX, e(2)-padY, e(3)+2*padX, e(4)+2*padY], ...
              'FaceColor','w','EdgeColor','none','Clipping','off');
    text(ax, x, y, str, 'FontSize', fs, 'FontName', fontname, ...
         'Color', color, 'Clipping','off');
end
