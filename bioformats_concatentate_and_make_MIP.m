
%
% ************************************************************************
% (c) 2023:
%       Miranda Hunter, White Lab, MSKCC
%           hunterm@mskcc.org | mirandavhunter@gmail.com
%
% Import .lif or .czi multi-channel timelapse confocal stacks, make maximum intensity projections, and save as multipage .tif.
% Optional: perform background subtraction using Gaussian blur method
% ************************************************************************

% Requires:
% Bio-Formats plugin
% saveastiff function from Matlab Central https://www.mathworks.com/matlabcentral/fileexchange/35684-multipage-tiff-stack


%% Setup

clear all
close all

data_dir = '/Volumes/whitelab/Lab Members/MirandaHunter/Microscopy/LSM880/230606_confinement_timelapse_siRNA/';

run_bg_subtraction = 0; % run background subtraction and save movie files?
if run_bg_subtraction
    gsig = 125; % sigma for Gaussian blur
end

%% Import data

% find movie files in folder
cd(data_dir);
movie_files = dir;
movie_file_names = {movie_files.name};
czis = movie_file_names(contains(movie_file_names, '.czi'));

for ii = 1:length(czis)

    % waitbar
    if ii == 1
        f = waitbar((ii-1)/length(czis), append('Processing movie ', string(ii), ' of ', string(length(czis)), '...'));
    else
        waitbar((ii-1)/length(czis), f, append('Processing movie ', string(ii), ' of ', string(length(czis)), '...'));
    end

    data = bfopen(append(data_dir, czis{ii}));

    % images within data are stored in this format:
    % each row corresponds to 1 series (containing all channels)
    % within each element of the first column (1,1 2,1 3,1 etc) there is an image file that corresponds to each individual image within the stack
    % images are stored as blue image (z-stack position 1), green image (z-stack position 1), red image (z-stack position 1), blue image (z-stack position 2) etc

    % Figure out parameters of data
    stack = data{1,1};

    n_timepoints = str2num(extractAfter(stack{1,2}, 'T=1/'));
    n_slices = extractBetween(stack{1,2}, 'Z=1/', ';');

    % dealing with situations where there is only 1 Z-slice (i.e FRAP)
    try
        n_slices = str2num(n_slices{1});
    catch
        fprintf('Only 1 Z-slice detected...\n')
        n_slices = 1;
    end

    [n_images_total,~] = size(data{1,1});
    n_channels = n_images_total ./ (n_timepoints*n_slices);

    im_names = stack(:,2);

    % Go time point by time point, create MIPs
    switch n_channels
        case 1

            movie_data = cell(n_timepoints,1);

            if run_bg_subtraction
                movie_bsub = cell(n_timepoints,1);
            end

            movie_channel1 = cell(n_timepoints,1);
            all_ims_ch1 = stack(1:n_timepoints,1);

            for jj = 1:n_timepoints

                timepoint_stack_ch1 = all_ims_ch1((n_slices*(jj-1)+1):(jj*n_slices),1);
                timepoint_stack_mat_ch1 = cat(3, timepoint_stack_ch1{:});
                mip_ch1 = max(timepoint_stack_mat_ch1, [], 3);
                movie_channel1{jj,1} = mip_ch1;

                clearvars('mip_ch*', 'timepoint_stack*');
            end

            % make 3D stack of all images per individual movie
            movie_stack_ch1 = cat(3, movie_channel1{:});

            % save
            cd(data_dir);
            clear options;
            options.overwrite = true;
            saveastiff(movie_stack_ch1, append(extractBefore(czis{ii}, '.czi'), '_MIP.tif'), options);

            clearvars('data')


        case 2
            movie_channel1 = cell(n_timepoints,1);
            movie_channel2 = cell(n_timepoints,1);

            if run_bg_subtraction
                movie_channel1_bsub = cell(n_timepoints,1);
                movie_channel2_bsub = cell(n_timepoints,1);
            end

            idx_ch1 = find(~cellfun(@isempty,(strfind(im_names, 'C=1/2'))));
            idx_ch2 = find(~cellfun(@isempty,(strfind(im_names, 'C=2/2'))));

            all_ims_ch1 = stack(idx_ch1,1);
            all_ims_ch2 = stack(idx_ch2,1);

            for jj = 1:n_timepoints

                timepoint_stack_ch1 = all_ims_ch1((n_slices*(jj-1)+1):(jj*n_slices),1);
                timepoint_stack_ch2 = all_ims_ch2((n_slices*(jj-1)+1):(jj*n_slices),1);

                % convert from cell to 3D matrix
                timepoint_stack_mat_ch1 = cat(3, timepoint_stack_ch1{:});
                timepoint_stack_mat_ch2 = cat(3, timepoint_stack_ch2{:});

                % make MIP
                mip_ch1 = max(timepoint_stack_mat_ch1, [], 3);
                mip_ch2 = max(timepoint_stack_mat_ch2, [], 3);

                movie_channel1{jj,1} = mip_ch1;
                movie_channel2{jj,1} = mip_ch2;

                % optional: background subtraction
                if run_bg_subtraction
                    ch1_gauss_im = imgaussfilt(mip_ch1, gsig);
                    ch2_gauss_im = imgaussfilt(mip_ch2, gsig);

                    movie_channel1_bsub{jj,1} = mip_ch1 - ch1_gauss_im;
                    movie_channel2_bsub{jj,1} = mip_ch2 - ch2_gauss_im;
                end

                clearvars('mip_ch*', 'timepoint_stack*');

            end

            % make 3D stack of all images per individual movie
            movie_stack_ch1 = cat(3, movie_channel1{:});
            movie_stack_ch2 = cat(3, movie_channel2{:});

            % save
            cd(data_dir);
            clear options;
            options.overwrite = true;
            saveastiff(movie_stack_ch1, append(extractBefore(czis{ii}, '.czi'), '_MIP_ch1.tif'), options);
            saveastiff(movie_stack_ch2, append(extractBefore(czis{ii}, '.czi'), '_MIP_ch2.tif'), options);

            if run_bg_subtraction
                movie_stack_ch1_bsub = cat(3, movie_channel1_bsub{:});
                movie_stack_ch2_bsub = cat(3, movie_channel2_bsub{:});

                cd(data_dir);
                saveastiff(movie_stack_ch1_bsub, append(extractBefore(czis{ii}, '.czi'), '_MIP_bsub_ch1.tif'), options);
                saveastiff(movie_stack_ch2_bsub, append(extractBefore(czis{ii}, '.czi'), '_MIP_bsub_ch2.tif'), options);
            end
            clearvars('data');
    end

end


waitbar(1, f, 'Done!');
pause(3)
close(f)

clear all




