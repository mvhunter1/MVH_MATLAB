% **********************************************************************************************************
% (c) 2022:
%       Miranda Hunter, White Lab, MSKCC
%           hunterm@mskcc.org | mirandavhunter@gmail.com
% 
% Script to analyze AFM results.
% - Import .tsv file with Young's modulus and height per pixel.
% - Convert to matrix to allow for plotting and segmentation as image.
% - Segment cell based on Young's modulus (since glass is much stiffer than the cell).
% - Segment nucleus based on cell height (since the nucleus is much taller than the rest of the cell.)
% - Calculate mean stiffness per pixel and plot.
% **********************************************************************************************************


clear all
close all

display_ims = 0;
display_thresh_ims = 0;

data_folder = '/Volumes/whitelab/Lab Members/MirandaHunter/AFM/221130_Ycompound/';

height_thresh = 12; % in µm
ym_thresh = 5; % in kPa


%% Import data

% find subfolders containing tsv file with results
cd(data_folder)
files = dir;
folders = files([files.isdir]);
folder_names = {folders.name};
folder_names = folder_names(~matches(folder_names, {'.', '..'}));

n_conditions = length(folder_names);
conditions_full = folder_names;
conditions = extractBefore(conditions_full, "_");

ym_all = {};
ym_all_thresh = {};
height_all = {};
xpos_all = {};
ypos_all = {};


% go through each folder and load data
for ii = 1:length(folder_names)
    folder_path = fullfile(data_folder, folder_names{ii});
    cd(folder_path)

    % again find each folder corresponding to each individual image/cell
    files_im = dir;
    folders_im = files_im([files_im.isdir]);
    folder_names_im = {folders_im.name};
    folder_names_im = folder_names_im(~matches(folder_names_im, {'.', '..'}));

    index = 0;
    % go into each folder and load TSV file
    for jj = 1:length(folder_names_im)

        fprintf(append('Loading data for image ', folder_names_im{jj}, ' from condition ', folder_names{ii}, '.\n'))

        image_folder_path = fullfile(folder_path, folder_names_im{jj});
        cd(image_folder_path);

        tsv_file = dir("*.tsv");
        tsv = tdfread(tsv_file.name);

        % extract relevant variables
        ym_Pa = tsv.Young0x27s_Modulus_0x5BPa0x5D;
        xpos = tsv.X_Position;
        ypos = tsv.Y_Position;
        height = tsv.Contact_Point_0x5Bm0x5D *1e6; % convert to µm


        % check if there are nans in the ym measurements. If so, throw out and move to next image.
        if any(isnan(ym_Pa))
            fprintf(append('NaNs detected in image ', folder_names_im{jj}, '. Skipping to next image.\n'))
            continue
        elseif (length(ypos) < 1024 || length(xpos) < 1024) % check if images are smaller than they're supposed to be.
            fprintf(append('Missing data detected in image ', folder_names_im{jj}, '. Skipping to next image.\n'))
            continue
        else
            index = index+1; % this preserves the structure of the cell array while removing empty elements.

            % for Youngs modulus, convert to kPa (original units are Pa) 
            ym_kPa = ym_Pa ./1000;

            % store data
            ym_all{ii,index} = ym_kPa;
            xpos_all{ii,index} = xpos;
            ypos_all{ii,index} = ypos;
            height_all{ii,index} = height;

        end
    end
end

[n_conditions, ~] = size(ym_all);
im_size = sqrt(length(ypos));


clearvars -except height_all n_conditions xpos_all ypos_all ym_all display_ims display_thresh_ims im_size conditions conditions_full height_thresh ym_thresh


%% Calculate number of samples per condition

n_samples = sum(~cellfun(@isempty, height_all),2)';

%% Make matrices corresponding to each image, so can plot as image instead of scatterplot.

fprintf('\nConverting data to matrices...\n')

% define X and Y coordinates to make matrix
xvals = repmat([1:im_size],1,im_size)';
yvals = flip(repelem([1:im_size],im_size)');

ym_matrices = {};
height_matrices = {};


for ii = 1:n_conditions
    for jj = 1:n_samples(ii)
        ym_zvals = ym_all{ii,jj};
        ym_matrices{ii,jj} = accumarray([xvals(:),yvals(:)],ym_zvals(:),[im_size im_size]);

        height_zvals = height_all{ii,jj};
        height_matrices{ii,jj} = accumarray([xvals(:),yvals(:)],height_zvals(:),[im_size im_size]);

        if display_ims
            figure; tiledlayout(1,2);
            nexttile; imagesc(ym_matrices{ii,jj}); colormap(gca, 'hot'); colorbar(gca); axis off
            nexttile; imagesc(height_matrices{ii,jj}); colormap(gca, 'parula'); colorbar(gca); axis off
        end
    end
end




%% Use YM to threshold height images.

fprintf('\nCreating cell masks...\n')

im_masks = {};
ym_all_thresh = {};
height_matrices_thresh = {};

for ii = 1:n_conditions
    for jj = 1:n_samples(ii)
        ym_data_plot = ym_matrices{ii,jj};

        ym_mask = ym_data_plot < ym_thresh;

        % use some morphological operations to clean up the thresholding.
        se = strel('square', 2);
        ym_mask = imerode(ym_mask, se);
        ym_mask = bwareafilt(ym_mask, [50 Inf]); 
        ym_mask = imdilate(ym_mask, se);

        height_matrices_thresh{ii,jj} = height_matrices{ii,jj} .* ym_mask;
        ym_matrices_thresh{ii,jj} = ym_matrices{ii,jj} .* ym_mask;
        im_masks{ii,jj} = ym_mask;

%         % add a slight erosion to remove edge (substrate) effect.
%         se2 = strel('square', 4);
%         ym_mask_erode = imerode(ym_mask, se2);
%         im_masks_erode{ii,jj} = ym_mask_erode;

        if display_thresh_ims
            figure; tiledlayout(1,3);
            nexttile; imagesc(ym_mask); colormap(gca, 'gray'); title('Mask'); axis off
            nexttile; imagesc(ym_matrices_thresh{ii,jj}); colormap(gca, 'hot'); colorbar(gca); title('Youngs modulus'); axis off; caxis(gca, [0 5]);
            nexttile; imagesc(height_matrices_thresh{ii,jj}); colormap(gca, 'parula'); colorbar(gca); title('Height'); axis off; caxis(gca, [0 20]);
        end
    end
end


%% Segment just the nucleus using height data

nuclear_thresh = {};
height_nuclear_thresh = {};
ym_nuclear_thresh = {};

for ii = 1:n_conditions
    for jj = 1:n_samples(ii)

        height_im = height_matrices_thresh{ii,jj};
        
        height_im_thresh = height_im > height_thresh;
        height_im_thresh = bwareafilt(height_im_thresh, [10 Inf]);
        %height_im_thresh = imclearborder(height_im_thresh);

        nuclear_thresh{ii,jj} = height_im_thresh;
        height_nuclear_thresh{ii,jj} = height_matrices{ii,jj} .* height_im_thresh;
        ym_nuclear_thresh{ii,jj} = ym_matrices{ii,jj} .* height_im_thresh;

        if display_thresh_ims
            figure; tiledlayout(1,2);
            nexttile; imagesc(height_im); colormap(gca, 'parula'); axis off; colorbar(gca); caxis([0 20]);
            nexttile; imagesc(height_im_thresh); colormap(gca, 'gray'); axis off;
        end
    end
end


%% Make figure displaying stiffness and nuclear mask for each condition

% cellular stiffness
for ii = 1:n_conditions
    figure;
    f = tiledlayout(4,7);
    for jj = 1:n_samples(ii)

        if jj < 8
            nexttile(jj); imagesc(im_masks{ii,jj}); colormap(gca, 'gray'); axis off
            nexttile(jj+7); imagesc(ym_matrices_thresh{ii,jj}); colormap(gca, 'hot'); axis off; caxis([0 5]);
        else
            nexttile(jj+7); imagesc(im_masks{ii,jj}); colormap(gca, 'gray'); axis off;
            nexttile(jj+14); imagesc(ym_matrices_thresh{ii,jj}); colormap(gca, 'hot'); axis off; caxis([0 5]);
        end
    end
    cb = colorbar('hot');
    cb.Layout.Tile = 'east';
    cb.FontSize = 20;
    cb.Color = 'white';
    title(f, conditions{ii}, 'FontWeight', 'normal');
end


% nuclear stiffness
for ii = 1:n_conditions
    figure;
    f = tiledlayout(4,7);
    for jj = 1:n_samples(ii)

        if jj < 8
            nexttile(jj); imagesc(nuclear_thresh{ii,jj}); colormap(gca, 'gray'); axis off
            nexttile(jj+7); imagesc(ym_nuclear_thresh{ii,jj}); colormap(gca, 'hot'); axis off; caxis([0 5]);
        else
            nexttile(jj+7); imagesc(nuclear_thresh{ii,jj}); colormap(gca, 'gray'); axis off;
            nexttile(jj+14); imagesc(ym_nuclear_thresh{ii,jj}); colormap(gca, 'hot'); axis off; caxis([0 5]);
        end
    end
    cb = colorbar('hot');
    cb.Layout.Tile = 'east';
    cb.FontSize = 20;
    cb.Color = 'white';
    title(f, conditions{ii}, 'FontWeight', 'normal');
end


%% Average stiffness values per cell.

fprintf('\n')
fprintf('Calculating stiffness per cell...\n')

mean_stiffness = nan(100, n_conditions);
mean_nuclear_stiffness = nan(100, n_conditions);

for ii = 1:n_conditions
    for jj = 1:n_samples(ii)
        
        ym_thresh_vals = ym_matrices_thresh{ii,jj};
        ym_thresh_vals = ym_thresh_vals(:);
        mean_stiffness(jj,ii) = mean(ym_thresh_vals(ym_thresh_vals > 0));

        ym_nuclear_vals = ym_nuclear_thresh{ii,jj};
        ym_nuclear_vals = ym_nuclear_vals(:);
        mean_nuclear_stiffness(jj,ii) = mean(ym_nuclear_vals(ym_nuclear_vals > 0));

    end
end


%% Plot 

cols = [172 172 172;
    191 0 53] ./256;


% do stats and add to title.
[~,p_cell] = ttest2(mean_stiffness(:,1), mean_stiffness(:,2));
[~,p_nuc] = ttest2(mean_nuclear_stiffness(:,1), mean_nuclear_stiffness(:,2));

figure; tiledlayout(1,2);

nexttile;
boxplotMVH(mean_stiffness, cols, conditions);
ylim([0 3]);
set(gca, 'YTick', (0:1:100));
ylabel({'Youngs modulus (kPa)'})
title('whole-cell stiffness');
subtitle(append('P = ', string(round(p_cell,5))));


% nuclear stiffness
nexttile;
boxplotMVH(mean_nuclear_stiffness, cols, conditions);
ylim([0 3]);
set(gca, 'YTick', (0:1:100));
ylabel({'Youngs modulus (kPa)'})
title('nuclear stiffness');
subtitle(append('P = ', string(round(p_nuc,5))));










