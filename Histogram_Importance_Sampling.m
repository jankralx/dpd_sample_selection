function [ us, used_hist ] = Histogram_Importance_Sampling(x, target_num_points, boundaries, method, reference)
% Histogram_Importance_Sampling function returns vector of selected points.
% The points are selected based on the histogram of points
% amplitudes. The points are selected such the output histogram would
% contains ideally target_num_points in all bins defined by given
% boundaries.
%
%   Inputs:
%   =======
%   x  - signal which is used for point selection based on abs(x)
%
%   target_num_points - target number of points in all bins of histogram
%                       with selected points
%
%   boundaries - boundaries of histogram bins
%
%   [method] - {'even', 'orig', 'ref', 'given', 'givex'} - goal for the output histogram. When
%   'even' is selected then the goal is to select points which create even
%   histogram with approx. same number of samples in all bins.
%
%   Returns:
%   ========
%
%   us - contains true if the corresponding sample is used, false
%   otherwise. It can be used directly for indexing variables.


% Authors: Jan Kral <kral.j@lit.cz>
% Date: 6.9.2018

if nargin < 4
    method = 'even';
end
if nargin < 6
    vicinity = 5;
end

[bin_cnts, ~, x_bins] = histcounts(x, boundaries);
bin_cnts = bin_cnts.';

if strcmp(method, 'orig')
    reference = x;
end

if strcmp(method, 'orig') || strcmp(method, 'ref')
    % output histogram should be similar to the histogram of the given
    % reference signal
    ref_cnts = histcounts(reference, boundaries);
    target_cnts = round(ref_cnts / (sum(ref_cnts)/target_num_points));
elseif strcmp(method, 'given') || strcmp(method, 'givex')
    % check that the reference histogram is compliant with given boundaries
    if length(boundaries) ~= length(reference) + 1
        error('Error: Size of the given histogram is different from given boundaries')
    end
    target_cnts_prec = reference / (sum(reference)/target_num_points);
    target_cnts = round(target_cnts_prec);
    
    % correction of target_cnts to have sum of target_num_points
    target_diff = sum(target_cnts) - target_num_points;
    if target_diff > 0
        % too much selected points
        target_cnts_prec_diff = target_cnts_prec - target_cnts;
        while target_diff ~= 0
            % add point to bin which has the most priority (is the closest
            % to the whole number)
            [~, ind] = min(target_cnts_prec_diff);
            target_cnts(ind) = target_cnts(ind) - 1;
            target_cnts_prec_diff(ind) = 0;
            target_diff = target_diff - 1;
        end
    elseif target_diff < 0
        % too less selected points
        target_cnts_prec_diff = target_cnts_prec - target_cnts;
        while target_diff ~= 0
            % add point to bin which has the most priority (is the closest
            % to the whole number)
            [~, ind] = max(target_cnts_prec_diff);
            target_cnts(ind) = target_cnts(ind) + 1;
            target_cnts_prec_diff(ind) = 0;
            target_diff = target_diff + 1;
        end
    end
    
    if sum(target_cnts) ~= target_num_points
        warning('Histogram acttarget_sumual number of points %i different from requested %i', sum(target_cnts), target_num_points);
    end
    
elseif strcmp(method, 'even')
    target_pts_in_bin = floor(target_num_points / length(bin_cnts));
    target_cnts = ones(size(bin_cnts))*target_pts_in_bin;
end

%% Check and modify number of target_cnts
% if number of sum of target counts differ from required

target_sum = sum(target_cnts);

if ~strcmp(method, 'givex')
    if target_sum > target_num_points
        % randomly select histogram bins whos values should be decreased
        sum_diff = abs(target_sum - target_num_points);
        rand_sel = randperm(length(target_cnts), sum_diff);
        target_cnts(rand_sel) = target_cnts(rand_sel) - 1;
    elseif target_sum < target_num_points
        % randomly select histogram bins whos values should be increased
        sum_diff = abs(target_sum - target_num_points);
        rand_sel = randperm(length(target_cnts), sum_diff);
        target_cnts(rand_sel) = target_cnts(rand_sel) + 1;
    end
else
    target_num_points = target_sum;
end

% check if the modification is not too drastic
if abs(target_sum - target_num_points) > 1
    warning('Modifying target histogram by more than 1 point: target_num_points = %d, actual num of points = %d', ...
        target_num_points, target_sum)
end

%% select points to the new histogram

dest_cnts = zeros(size(bin_cnts));
us = zeros(target_num_points, 1);
us_len = 0;
% go through all bins and randomly select points into the new bin
for i = 1:length(bin_cnts)
    if bin_cnts(i) <= target_cnts(i)
        % there is too few points in the bin, take them all
        select = find(x_bins == i);
    else
        % enough points, we need to select only few of them
        pts_indx = find(x_bins == i);
        select = pts_indx(randperm(length(pts_indx),target_cnts(i)));
    end
    us(us_len+1:us_len+length(select)) = select;
    dest_cnts(i) = length(select);
    us_len = us_len + dest_cnts(i);
end

%% Correct the number of points in the output histogram

% go throught output histogram and add samples if there is insuficient
% number of points in certain bins
for i = 1:length(bin_cnts)
    if dest_cnts(i) < target_cnts(i)
        % not enough points in the target bin, try to find new points in
        % neighbour bins
        added_pts = 0;
        ctn_diff = target_cnts(i) - dest_cnts(i);
        k = 1;
        while added_pts < ctn_diff && k < 2*length(bin_cnts)
            k = k + 1;
            d = floor(k/2);
            sign = -1*mod(k,2);
            bin_ind = i + sign*d;
            if bin_ind < 1 || bin_ind > length(bin_cnts)
                continue;
            end
            % add missing number of points from bin_ind to output histogram
            pts_indx = find(x_bins == bin_ind);     % selects points which belongs to the bin_ind bin
            pts_indx = setdiff(pts_indx, us);       % from the selected points remove those already in the output
            if length(pts_indx) <= ctn_diff - added_pts
                % not enough points in the newly selected bin
                % add all of them
                select = pts_indx;
            else
                select = pts_indx(randperm(length(pts_indx),ctn_diff-added_pts));    % randomly select the new points
            end
            us(us_len+1:us_len+length(select)) = select;
            us_len = us_len + length(select);
            added_pts = added_pts + length(select);
        end
    end
end
us = us(1:us_len);

% sort undersampling matrix and check its properties
us = sort(us);

% check the length
if length(us) ~= target_num_points
%    warning('WARNING: Output of histogram method contains %i different from requested %i', length(us), target_num_points);
end

% % check that the selected elements are not close to each other
% if min(diff(us)) < 5
%     warning('WARNING: Output of histogram method contains points close to each other');
% end

% save used histogram to be returned
used_hist = target_cnts;

% % plot histograms
% figure(21);
% histogram(x, boundaries);
% hold on;
% histogram(x(us), boundaries);
% hold off;
% axis([boundaries(1) boundaries(end), 0 ceil(max(target_cnts)*1.4)]);

return

%% this is old unused method
%  The mothod did not provide the required results as the output histogram
%  was only approximation of the required one.
[bin_cnts, ~, x_bins] = histcounts(abs(x), boundaries);
bin_cnts = bin_cnts.';
target_pts_in_bin = target_num_points / (length(boundaries)-1);
bin_probability = target_pts_in_bin ./ bin_cnts;
x_prob = bin_probability(x_bins);
randv = rand(size(x));
us = randv < x_prob;

% plot histograms
figure(21);
histogram(abs(x), boundaries);
hold on;
histogram(abs(x(us)), boundaries);
hold off;
axis([boundaries(1) boundaries(end), 0 ceil(target_pts_in_bin*1.4)]);

end

