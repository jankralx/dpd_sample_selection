close all;
addpath('mutils');

% define default colors  for plots
colors = [
         0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0         0.5       0
    0.4940    0.1840    0.5560
    0.9290    0.6940    0.1250
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840
    ];

markers = ['x'; 'o'; 's'; '^'; '+'; '*'];
line_styles = {'-','--',':','-.'};
no_dpd_line_width = 1.5;

is_no_dpd = 0;
is_annotation = 0;
is_pdf_export = 0;
file_name = 'results_01';

% check if nmse_all exists in variables, if not, load it from results file
if ~exist('nmse_all', 'var')
    % check if there exists results_cache file, if exists, load it
    if exist(sprintf('%s_cache.mat', file_name), 'file')
	load(sprintf('%s_cache.mat', file_name));
    else
	load(sprintf('%s.mat', file_name));

	% extract needed data from data saved by calculation
	evm_all = zeros(length(num_samples_simrange), pars.num_arch, pars.num_prob_rep, pars.num_iter);
	nmse_all = zeros(length(num_samples_simrange), pars.num_arch, pars.num_prob_rep, pars.num_iter);
	acpr_all = zeros(length(num_samples_simrange), pars.num_arch, pars.num_prob_rep, pars.num_iter, 2, 2);

	for i = 1:length(num_samples_simrange)
	    for m = 1:pars.num_prob_rep
            for k = 1:pars.num_arch
                if isempty(arch_res{i,m})
                    continue
                end
                % EVM
                evm_all(i,k,m,:) = arch_res{i,m}(k).evm;
                % NMSE
                nmse_all(i,k,m,:) = [arch_res{i,m}(k).res(:).nmse].';
                % ACPR
                acpr_all(i,k,m,:,:,:) = arch_res{i,m}(k).acpr(:,:,:);
            end
	    end
	end

	% save extracted data into cache file
	save(sprintf('%s_cache.mat', file_name), 'evm_all', 'nmse_all', 'acpr_all', 'arch', 'pars', 'num_samples_simrange', '-v7.3');
    end
end

if is_no_dpd
    % load results for PA without DPD and for signal itself
    load(sprintf('%s_nodpd.mat', file_name), 'ares_no_dpd', 'ares_no_pa');
    % calculate metrics for PA without dpd
    nmse_no_dpd = Avg_dB([ares_no_dpd.res(:).nmse], 2, 10);
    evm_no_dpd = rms(ares_no_dpd.evm);
    evm_no_pa = rms(ares_no_pa.evm);
    acpr_no_dpd = reshape(Avg_dB(ares_no_dpd.acpr, 1, 10), 2, 2);
end

% convergence interval definitions
ci_percentage = 95;
ci_len = round((pars.num_prob_rep - pars.num_prob_rep * ci_percentage / 100)/2);

% calculate averages in avg_fft_after interval
evm_avg = zeros(length(num_samples_simrange), pars.num_prob_rep, pars.num_arch);
nmse_avg = zeros(length(num_samples_simrange), pars.num_prob_rep, pars.num_arch);
acpr_avg = zeros(length(num_samples_simrange), pars.num_prob_rep, pars.num_arch, 4);

for i = 1:length(num_samples_simrange)
    for m = 1:pars.num_prob_rep
        for k = 1:pars.num_arch
            % EVM
            evm_avg(i,m,k) = rms(evm_all(i,k,m,pars.res.avg_fft_after_iter:end));
            % NMSE
            nmserr = nmse_all(i,k,m,:);
            nmserr = nmserr(:);
            nmse_avg(i,m,k) = Avg_dB(nmserr(pars.res.avg_fft_after_iter:end), 1, 10);
            % ACPR
            acprr = Avg_dB(acpr_all(i,k,m,pars.res.avg_fft_after_iter:end,:,:),4,10);
            acprr = reshape(acprr,2,2);
            acpr_avg(i,m,k,:) = reshape(acprr.',4,1);
        end
    end
end

% detect runs that did not converge
converged = nmse_avg < 0;


%% NMSE
fig = figure(1);
pos = get(gcf, 'Position');
pos(3:4) = [500 375];
set(gcf, 'Position', pos);
%leg_txt = {'Sample Block', 'Decimation 1/F', 'EDH', 'GOH', 'QR', 'GS', 'GS+QR'};
leg_txt = {'Conv. DPD', 'US', 'EDH', 'GOH', 'QRS', 'GSS'};

if is_no_dpd
    % plot a line for no dpd value
    plot([10; 100e3], repmat(nmse_no_dpd,2,1), 'k--', 'LineWidth', no_dpd_line_width);
end

nmse_avg2 = zeros(length(num_samples_simrange),1);
nmse_ci = zeros(2,length(num_samples_simrange));
for k = 1:pars.num_arch-1
    for i = 1:length(num_samples_simrange)
        % calculate average of all probability runs which has converged
        nmse_avg2(i) = Avg_dB(nmse_avg(i,converged(i,:,k),k), 2, 10);
        try
            nmserr = sort(nmse_avg(i,converged(i,:,k),k), 2);
            nmserr = nmserr(1+ci_len:end-ci_len);
            nmse_ci(1,i) = min(nmserr);
            nmse_ci(2,i) = max(nmserr);
        catch
            nmse_ci(:,i) = NaN;
            nmse_avg2(i) = NaN;
        end
    end
    nmse_range = ~isnan(nmse_avg2);
    nmse_plot = [nmse_avg2, nmse_ci(1,:).', nmse_ci(2,:).'];
    hold on;
    plot_ci(num_samples_simrange(nmse_range), nmse_plot(nmse_range,:), ...
        'PatchColor', colors(k,:), 'PatchAlpha', 0.1, ...
        'MainLineWidth', 2, 'MainLineStyle', '-', 'MainLineColor', colors(k,:), ...
        'LineWidth', 1.5, 'LineStyle','--', 'LineColor', colors(k,:));
    hold off;

end
set(gca, 'XScale', 'log');
ylim([-36 -10]);
grid on;
if is_annotation
    % draw annotation arrows with text into figure
    arrow_normlen = [1/20 0 1/15 0];
    pos = [xy2norm('x', 700, 700) xy2norm('y', -23.09, -23.09)] + arrow_normlen;
    annotation('textarrow', pos(1:2), pos(3:4), 'String', leg_txt{1});
    pos = [xy2norm('x', 300, 300) xy2norm('y', -22.03, -22.03)] + 2.1*arrow_normlen;
    annotation('textarrow', pos(1:2), pos(3:4), 'String', leg_txt{2});
    pos = [xy2norm('x', 50, 50) xy2norm('y', -34.13, -34.13)] + 0.5*arrow_normlen;
    annotation('textarrow', pos(1:2), pos(3:4), 'String', leg_txt{3});
    pos = [xy2norm('x', 100, 100) xy2norm('y', -36.3, -36.3)] + arrow_normlen;
    annotation('textarrow', pos(1:2), pos(3:4), 'String', leg_txt{4});
end

% disable unwanted lines of legend
h2 = get(gca, 'Children');
for i = 1:length(h2)
    if mod(i+1,4)
         h2(i).HandleVisibility = 'off';
    end
end
h2(end).HandleVisibility = 'on';
legend([{'No DPD'} leg_txt]);
xlabel('Number of Selected Samples $N$ (-)');
ylabel('NMSE (dB)');
ApplyFigureSettings(fig);
if is_pdf_export
    saveas(gcf, 'figures/hondri_nmse.pdf');
    saveas(gcf, 'figures/hondri_nmse.fig');
end
%% EVM
fig = figure(2);
pos = get(gcf, 'Position');
pos(3:4) = [500 375];
set(gcf, 'Position', pos);

if is_no_dpd
    % plot a line for no pa value
%     plot([10; 100e3], repmat(evm_no_pa*100,2,1), '--', 'Color', [0.5 0.5 0.5]);
    % plot a line for no dpd value
%     hold on;
    plot([10; 100e3], repmat(evm_no_dpd*100,2,1), 'k--', 'LineWidth', no_dpd_line_width);
%     hold off;
end

evm_avg2 = zeros(length(num_samples_simrange),1);
evm_ci = zeros(2,length(num_samples_simrange));
for k = 1:pars.num_arch-1
    for i = 1:length(num_samples_simrange)
        % calculate average of all probability runs which has converged
        evm_avg2(i) = rms(evm_avg(i,converged(i,:,k),k), 2);
        try
            evmr = sort(evm_avg(i,converged(i,:,k),k), 2);
            evmr = evmr(1+ci_len:end-ci_len);
            evm_ci(1,i) = min(evmr);
            evm_ci(2,i) = max(evmr);
        catch
            evm_ci(:,i) = NaN;
            evm_avg2(i) = NaN;
        end
    end
    evm_range = ~isnan(evm_avg2);
    evm_plot = [evm_avg2, evm_ci(1,:).', evm_ci(2,:).'];

    hold on;
    plot_ci(num_samples_simrange(evm_range), 100*evm_plot(evm_range,:), ...
        'PatchColor', colors(k,:), 'PatchAlpha', 0.1, ...
        'MainLineWidth', 2, 'MainLineStyle', '-', 'MainLineColor', colors(k,:), ...
        'LineWidth', 1.5, 'LineStyle','--', 'LineColor', colors(k,:));
    hold off;
end
set(gca, 'XScale', 'log');
ylim([1.3 6]);
grid on;
if is_annotation
    % draw annotation arrows with text into figure
    arrow_normlen = [1/20 0 1/15 0];
    pos = [xy2norm('x', 2000, 2000) xy2norm('y', 4.593, 4.593)] + arrow_normlen;
    annotation('textarrow', pos(1:2), pos(3:4), 'String', leg_txt{1});
    pos = [xy2norm('x', 1000, 1000) xy2norm('y', 4.586, 4.586)] + 2.2*arrow_normlen;
    annotation('textarrow', pos(1:2), pos(3:4), 'String', leg_txt{2});
    pos = [xy2norm('x', 50, 50) xy2norm('y', 4.152, 4.152)] + arrow_normlen;
    annotation('textarrow', pos(1:2), pos(3:4), 'String', leg_txt{3});
    pos = [xy2norm('x', 100, 100) xy2norm('y', 4.059, 4.059)] + arrow_normlen;
    annotation('textarrow', pos(1:2), pos(3:4), 'String', leg_txt{4});
end

% disable unwanted lines of legend
h2 = get(gca, 'Children');
for i = 1:length(h2)
    if mod(i+1,4)
        h2(i).HandleVisibility = 'off';
    end
end
h2(end).HandleVisibility = 'on';
legend([{'No DPD'} leg_txt]);
xlabel('Number of Selected Samples $N$ (-)');
ylabel('EVM (\%)');
ApplyFigureSettings(fig);
if is_pdf_export
    saveas(gcf, 'figures/hondri_evm.pdf');
    saveas(gcf, 'figures/hondri_evm.fig');
end

% % probability of convergence
% figure(3);
% for k = 1:pars.num_arch
%     hold on;
%     plot(num_samples_simrange, sum(converged(:,:,k),2), ...
%         'LineStyle', line_styles{1}, 'Color', colors(k,:), ...
%         'DisplayName', sprintf('Avg-%s',arch(k).name));
%     hold off;
% end
% set(gca, 'XScale', 'log');
% title('Convergence probability');
% xlabel('Number of Selected Samples $N$ (-)');
% ylabel('Probability (%)');
% legend('show');


%% ACPR - 1st adjacent channel
fig = figure(3);
pos = get(gcf, 'Position');
pos(3:4) = [500 375];
set(gcf, 'Position', pos);

if is_no_dpd
    % plot a line for no dpd value
    acpr_no_dpd2 = Avg_dB(acpr_no_dpd, 2, 10);
    plot([10; 100e3], repmat(acpr_no_dpd2(1),2,1), 'k--', 'LineWidth', no_dpd_line_width);
end

acpr_avg2 = zeros(length(num_samples_simrange),1);
acpr_ci = zeros(2,length(num_samples_simrange));
for k = 1:pars.num_arch-1
    for i = 1:length(num_samples_simrange)
        % calculate average of all probability runs which has converged
        try
            acprr = acpr_avg(i,converged(i,:,k),k,1:2);
            acprr = reshape(acprr, numel(acprr), 1);
            acprr = sort(acprr);
            acprr = acprr(1+2*ci_len:end-2*ci_len);

            acpr_avg2(i) = Avg_dB(acprr, 1, 10);
            acpr_ci(1,i) = min(acprr);
            acpr_ci(2,i) = max(acprr);
        catch
            acpr_ci(:,i) = NaN;
            acpr_avg2(i) = NaN;
        end
    end
    acpr_range = ~isnan(acpr_avg2);
    acpr_plot = [acpr_avg2, acpr_ci(1,:).', acpr_ci(2,:).'];

    hold on;
    plot_ci(num_samples_simrange(acpr_range), acpr_plot(acpr_range,:), ...
        'PatchColor', colors(k,:), 'PatchAlpha', 0.1, ...
        'MainLineWidth', 2, 'MainLineStyle', '-', 'MainLineColor', colors(k,:), ...
        'LineWidth', 1.5, 'LineStyle','--', 'LineColor', colors(k,:));
    hold off;
end
set(gca, 'XScale', 'log');
ylim([-45 -15]);
grid on;

if is_annotation
    % draw annotation arrows with text into figure
    arrow_normlen = [1/20 0 1/15 0];
    pos = [xy2norm('x', 1000, 1000) xy2norm('y', -33.14, -33.14)] + arrow_normlen;
    annotation('textarrow', pos(1:2), pos(3:4), 'String', leg_txt{1});
    pos = [xy2norm('x', 400, 400) xy2norm('y', -31.03, -31.03)] + 2*arrow_normlen;
    annotation('textarrow', pos(1:2), pos(3:4), 'String', leg_txt{2});
    pos = [xy2norm('x', 50, 50) xy2norm('y', -42.25, -42.25)] + arrow_normlen;
    annotation('textarrow', pos(1:2), pos(3:4), 'String', leg_txt{3});
    pos = [xy2norm('x', 100, 100) xy2norm('y', -44.03, -44.03)] + 1.3*arrow_normlen;
    annotation('textarrow', pos(1:2), pos(3:4), 'String', leg_txt{4});
end

% disable unwanted lines of legend
h2 = get(gca, 'Children');
for i = 1:length(h2)
    if mod(i+1,4)
        h2(i).HandleVisibility = 'off';
    end
end
h2(end).HandleVisibility = 'on';
legend([{'No DPD'} leg_txt]);

% title('ACPR in 1st Adjacent Channel');
xlabel('Number of Selected Samples $N$ (-)');
ylabel('ACPR (dB)');
ApplyFigureSettings(fig);
if is_pdf_export
    saveas(gcf, 'figures/hondri_acpr1.pdf');
    saveas(gcf, 'figures/hondri_acpr1.fig');
end


%% ACPR - 2nd adjacent channel
fig = figure(4);
pos = get(gcf, 'Position');
pos(3:4) = [500 375];
set(gcf, 'Position', pos);

if is_no_dpd
    % plot a line for no dpd value
    plot([10; 100e3], repmat(acpr_no_dpd2(2),2,1), 'k--', 'LineWidth', no_dpd_line_width);
end

acpr_avg2 = zeros(length(num_samples_simrange),1);
acpr_ci = zeros(2,length(num_samples_simrange));
for k = 1:pars.num_arch
    for i = 1:length(num_samples_simrange)
        % calculate average of all probability runs which has converged
        try
            acprr = acpr_avg(i,converged(i,:,k),k,3:4);
            acprr = reshape(acprr, numel(acprr), 1);
            acprr = sort(acprr);
            acprr = acprr(1+2*ci_len:end-2*ci_len);

            acpr_avg2(i) = Avg_dB(acprr, 1, 10);
            acpr_ci(1,i) = min(acprr);
            acpr_ci(2,i) = max(acprr);
        catch
            acpr_ci(:,i) = NaN;
            acpr_avg2(i) = NaN;
        end
    end
    acpr_range = ~isnan(acpr_avg2);
    acpr_plot = [acpr_avg2, acpr_ci(1,:).', acpr_ci(2,:).'];
%     acpr_plot = [acpr_avg2, acpr_avg2, acpr_avg2];

    hold on;
    plot_ci(num_samples_simrange(acpr_range), acpr_plot(acpr_range,:), ...
        'PatchColor', colors(k,:), 'PatchAlpha', 0.1, ...
        'MainLineWidth', 2, 'MainLineStyle', '-', 'MainLineColor', colors(k,:), ...
        'LineWidth', 1.5, 'LineStyle','--', 'LineColor', colors(k,:));
    hold off;
end
set(gca, 'XScale', 'log');
ylim([-60 -10]);
grid on;
if is_annotation
    % draw annotation arrows with text into figure
    arrow_normlen = [1/20 0 1/15 0];
    pos = [xy2norm('x', 600, 600) xy2norm('y', -40.61, -40.61)] + 0.8*arrow_normlen;
    annotation('textarrow', pos(1:2), pos(3:4), 'String', leg_txt{1});
    pos = [xy2norm('x', 310, 310) xy2norm('y', -40.64, -40.64)] + 2*arrow_normlen;
    annotation('textarrow', pos(1:2), pos(3:4), 'String', leg_txt{2});
    arrow_normlen = [1/20 0 1/10 0];
    pos = [xy2norm('x', 40, 40) xy2norm('y', -52.06, -52.06)] + 1.1*arrow_normlen;
    h = annotation('textarrow', pos(1:2), pos(3:4), 'String', leg_txt{3});
    pos = [xy2norm('x', 80, 80) xy2norm('y', -51.34, -51.34)] + 0.5*arrow_normlen;
    annotation('textarrow', pos(1:2), pos(3:4), 'String', leg_txt{4});
end

% disable unwanted lines of legend
h2 = get(gca, 'Children');
for i = 1:length(h2)
    if mod(i+1,4)
        h2(i).HandleVisibility = 'off';
    end
end
h2(end).HandleVisibility = 'on';
legend([{'No DPD'} leg_txt]);

% title('ACPR in 1st Adjacent Channel');
xlabel('Number of Selected Samples $N$ (-)');
ylabel('ACPR (dB)');
ApplyFigureSettings(fig);
if is_pdf_export
    saveas(gcf, 'figures/hondri_acpr2.pdf');
    saveas(gcf, 'figures/hondri_acpr2.fig');
end

% figure(9);
% i = length(num_samples_simrange); % find(num_samples_simrange == 300);
% k = 2;
% for m = 1:pars.num_prob_rep
%     % NMSE
%     evm = arch_res{i,m}(k).evm;
%     hold on;
%     plot(evm, ...
%         'LineStyle', line_styles{1}, 'Color', colors(mod(m-1,size(colors,2))+1,:), ...
%         'DisplayName', sprintf('%s',arch(k).name));
%     hold off;
% end

%% Calculate and plot histograms for EVM, NMSE, ACPR

arch_names = {'Block', '1/N', 'Hist-Even', 'Hist-GA'};

% InteractiveHist(num_samples_simrange, arch_names, 100*evm_avg, 20);
% InteractiveHist(num_samples_simrange, arch_names, nmse_avg, 20);


%%

arch_names = {'Sample Block', 'Decimation 1/N', 'EDH', 'GOH', 'QR', 'GS', 'GS+QR'};

% InteractivePlot(num_samples_simrange, arch_names, 100*evm_all);
% InteractivePlot(num_samples_simrange, arch_names, nmse_all);


function InteractivePlot(pars1, pars2, y)
    % create a figure
    f = figure;
    ax = axes('Parent', f, 'position', [0.13 0.39  0.77 0.54]);
    for rep = 1:size(y,3)
        y_sel = y(1,1,rep,:);
        y_sel = reshape(y_sel, numel(y_sel), 1);
        hold on;
        h(rep) = plot(ax, y_sel);
        hold off;
    end

    bgcolor = f.Color;

    % add slider to control selection of data which is plot
    a = uicontrol('Parent',f,'Style','slider','Units','normal','Position',[0.17,0.25,0.7,0.05],...
                  'value', 1, 'min', 1, 'max', length(pars1));

    as = uicontrol('Parent',f,'Style','text','Units','normal','Position',[0,0.2,1,0.05],...
                   'String',sprintf('Num Samp = %d', min(pars1)),'BackgroundColor',bgcolor);

    b = uicontrol('Parent',f,'Style','slider','Units','normal','Position',[0.17,0.1,0.7,0.05],...
                  'value', 1, 'min', 1, 'max', length(pars2));

    bs = uicontrol('Parent',f,'Style','text','Units','normal','Position',[0,0.05,1,0.05],...
                   'String',sprintf('Arch Name = %s', pars2{1}),'BackgroundColor',bgcolor);

    a.Callback = @(es,ed) UpdatePlotData(h, pars1, pars2, y, a, as, b, bs);
    b.Callback = @(es,ed) UpdatePlotData(h, pars1, pars2, y, a, as, b, bs);
end

function ind = UpdatePlotData(h, pars1, pars2, y, sel1, sels1, sel2, sels2)
    ind = round(sel1.Value);
    ind2 = round(sel2.Value);
    for rep = 1:size(y,3)
        y_sel = y(ind,ind2,rep,:);
        y_sel = reshape(y_sel, numel(y_sel), 1);
        set(h(rep), 'YData', y_sel);
    end
    set(sels1, 'String', sprintf('Num Samp = %d', pars1(ind)));
    set(sels2, 'String', sprintf('Arch Name = %s', pars2{ind2}));
end



function InteractiveHist(pars1, pars2, y, x)
    % create a figure
    f = figure;
    ax = axes('Parent', f, 'position', [0.13 0.39  0.77 0.54]);
    y_sel = y(1,:,1);
    y_sel = reshape(y_sel, numel(y_sel), 1);
    h = histogram(ax, y_sel, x);
    %h(rep) = bar(ax, y_sel, x);

    bgcolor = f.Color;

    % add slider to control selection of data which is plot
    a = uicontrol('Parent',f,'Style','slider','Units','normal','Position',[0.17,0.25,0.7,0.05],...
                  'value', 1, 'min', 1, 'max', length(pars1));

    as = uicontrol('Parent',f,'Style','text','Units','normal','Position',[0,0.2,1,0.05],...
                   'String',sprintf('Num Samp = %d', min(pars1)),'BackgroundColor',bgcolor);

    b = uicontrol('Parent',f,'Style','slider','Units','normal','Position',[0.17,0.1,0.7,0.05],...
                  'value', 1, 'min', 1, 'max', length(pars2));

    bs = uicontrol('Parent',f,'Style','text','Units','normal','Position',[0,0.05,1,0.05],...
                   'String',sprintf('Arch Name = %s', pars2{1}),'BackgroundColor',bgcolor);

    a.Callback = @(es,ed) UpdateHistData(ax, pars1, pars2, y, x, a, as, b, bs);
    b.Callback = @(es,ed) UpdateHistData(ax, pars1, pars2, y, x, a, as, b, bs);
end

function ind = UpdateHistData(ax, pars1, pars2, y, x, sel1, sels1, sel2, sels2)
    ind = round(sel1.Value);
    ind2 = round(sel2.Value);

    y_sel = y(ind,:,ind2);
    y_sel = reshape(y_sel, numel(y_sel), 1);
%     cnts = histcounts(y_sel, x);
%     set(h, 'YData', cnts);
    histogram(ax, y_sel, x);

    set(sels1, 'String', sprintf('Num Samp = %d', pars1(ind)));
    set(sels2, 'String', sprintf('Arch Name = %s', pars2{ind2}));
end
