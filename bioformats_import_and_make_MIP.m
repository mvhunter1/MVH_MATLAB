% 
% ************************************************************************
% (c) 2022:
%       Miranda Hunter, White Lab, MSKCC
%           hunterm@mskcc.org | mirandavhunter@gmail.com
% 
% Import .lif or .czi multi-channel image stacks, make maximum intensity projections, and save as .tif.
% ************************************************************************

% Requires: Bio-Formats plugin

%% Setup

clear all
close all

% if multiple images are saved in a folder (leave empty otherwise):
path_to_folder = {};

% if multiple images are stored within 1 file (leave empty otherwise):
path_to_file = {'/Volumes/whitelab/Lab Members/MirandaHunter/Microscopy/SP5/220411_hmgb2/CRISPR_12wk_1_slide1_tumor488_hmgb2555_40X_zoom.lif'};

save_MIPs = 1; % save MIPs as tif?
n_channels = 3; % number of channels acquired


%% Check if filetype is lif or czi

if ~isempty(path_to_file)
    if endsWith(path_to_file{1}, '.lif')
        filetype = '.lif';
    elseif endsWith(path_to_file{1}, '.czi')
        filetype = '.czi';
    else 
        error('Filetype must be .lif or .czi.')
    end
elseif ~isempty(path_to_folder)
    cd(path_to_folder{1})
    folder_files = dir;
    folder_files = {folder_files.name};
    folder_files = folder_files(endsWith(folder_files, {'.czi', '.lif'}));
    if unique(endsWith(folder_files, '.czi')) 
        filetype = '.czi';
    elseif unique(endsWith(folder_files, '.lif'))
        filetype = '.lif';
    else
        error('Filetype must be .lif or .czi.')
    end
end

if filetype == '.lif'
    fprintf('.lif file type detected.\n')
elseif filetype == '.czi'
    fprintf('.czi file type detected.\n')
end


%% Import file and save each series

ims_405 = {};
ims_488 = {};
ims_555 = {};
if n_channels == 4
    ims_647 = {};
end

show_ims = 0;

% for multiple .czi files stored in 1 folder
if ~isempty(path_to_folder) 

    % find image files within folder
    im_folder = path_to_folder{1};
    im_paths = fullfile(im_folder, folder_files);
    
    % determine name of experiment (i.e name of folder)
    idx = regexp(im_folder, '/');
    if endsWith(im_folder, '/')
        idx_start = idx(length(idx)-1);
        experiment_name = extractAfter(im_folder, idx_start);
        experiment_name = experiment_name(1:end-1);
    else 
        idx_start = idx(end);
        experiment_name = extractAfter(im_folder, idx_start);
    end
    
    for ii = 1:length(im_paths)
        
        % waitbar
        if ii == 1
            f = waitbar((ii-1)/length(im_paths), append('Processing image ', string(ii), ' of ', string(length(im_paths)), '...'));
        else
            waitbar((ii-1)/length(im_paths), f, append('Processing image ', string(ii), ' of ', string(length(im_paths)), '...'));
        end
        
        file_import = im_paths{ii};
        data = bfopen(file_import);
        
        % images within data are stored in this format:
        % each row corresponds to 1 series (containing all channels)
        % within each element of the first column (1,1 2,1 3,1 etc) there is an image file that corresponds to each individual image within the stack
        % images are stored as blue image (z-stack position 1), green image (z-stack position 1), red image (z-stack position 1), blue image (z-stack position 2) etc
        
        % determine number of images within the series
        [n_series,~] = size(data);
        
        for jj = 1:n_series
            
            series = data{jj,1};
            series_names = series(:,2);
            
            % find the images corresponding to each channel
            if n_channels == 3
                idx_405 = find(~cellfun(@isempty,(strfind(series_names, 'C=1/3'))));
                idx_488 = find(~cellfun(@isempty,(strfind(series_names, 'C=2/3'))));
                idx_555 = find(~cellfun(@isempty,(strfind(series_names, 'C=3/3'))));
            elseif n_channels == 4
                idx_405 = find(~cellfun(@isempty,(strfind(series_names, 'C=1/4'))));
                idx_488 = find(~cellfun(@isempty,(strfind(series_names, 'C=2/4'))));
                idx_555 = find(~cellfun(@isempty,(strfind(series_names, 'C=3/4'))));
                idx_647 = find(~cellfun(@isempty,(strfind(series_names, 'C=4/4'))));
            end
            if isempty(idx_647)
                error('Selected number of channels not detected. Check to make sure n_channels is correct.')
            end
            
            % move the images for each channel into their own stack
            im_stack_405 = series(idx_405,1);
            im_stack_488 = series(idx_488,1);
            im_stack_555 = series(idx_555,1);
            
            % convert from cell to 3D matrix
            im_stack_405 = cat(3,im_stack_405{:});
            im_stack_488 = cat(3,im_stack_488{:});
            im_stack_555 = cat(3,im_stack_555{:});
            
            % make MIP
            mip_405 = max(im_stack_405, [], 3);
            mip_488 = max(im_stack_488, [], 3);
            mip_555 = max(im_stack_555, [], 3);
            
            if n_channels == 4
                im_stack_647 = series(idx_647,1);
                im_stack_647 = cat(3, im_stack_647{:});
                mip_647 = max(im_stack_647, [], 3);
            end
            
            if show_ims
                figure;
                subplot(1,3,1); imagesc(mip_405); axis off
                subplot(1,3,2); imagesc(mip_488); axis off
                subplot(1,3,3); imagesc(mip_555); axis off
            end
            
            if n_channels == 3 || n_channels == 4
                % store images in MATLAB
                index = length(ims_405) + 1;
                ims_405{index} = mip_405;
                ims_488{index} = mip_488;
                ims_555{index} = mip_555;
                if n_channels == 4
                    ims_647{index} = mip_647;
                end
            end
            
%              % extract metadata and create metadata object for save
%             %metadata_im = createMinimalOMEXMLMetadata(mip_405);
%             metadata = data{1,4};
%             pixelSizeX = metadata.getPixelsPhysicalSizeX(0).value();
%             pixelSizeX = pixelSizeX.doubleValue();
%             pixelSizeY = metadata.getPixelsPhysicalSizeY(0).value();
%             pixelSizeY = pixelSizeY.doubleValue();
%             
%             %pixelSize = ome.units.quantity.Length(java.lang.Double(pixelSizeX), ome.units.UNITS.MICROMETER);
%             %metadata_im.setPixelsPhysicalSizeX(pixelSize, 0);
%             %metadata_im.setPixelsPhysicalSizeY(pixelSize, 0);
%             
%             [px_x, px_y] = size(mip_405);
%             length_x = px_x * pixelSizeX;
%             length_y = px_y * pixelSizeY;
%             
            % write MIP to folder
            cd(im_folder)
            im_name = extractAfter(file_import, append(experiment_name, '/'));
            im_name = extractBefore(im_name, filetype);
         
            imwrite(mip_405, append('max_', im_name, '_405.tiff'))
            imwrite(mip_488, append('max_', im_name, '_488.tiff'))
            imwrite(mip_555, append('max_', im_name, '_555.tiff'))
            if n_channels == 4
                imwrite(mip_647, append('max_', im_name, '_647.tiff'));
            end
            
        end
    end
     waitbar(1, f, 'Done!');
     pause(3)
     close(f)

% for multiple images stored within a single file
elseif ~isempty(path_to_file)
    
    data = bfopen(path_to_file{1});
    
    % images within data are stored in this format:
    % each row corresponds to 1 series (containing all channels)
    % within each element of the first column (1,1 2,1 3,1 etc) there is an image file that corresponds to each individual image within the stack
    % images are stored as blue image (z-stack position 1), green image (z-stack position 1), red image (z-stack position 1), blue image (z-stack position 2) etc
    
    % determine number of images within the series
    [n_series,~] = size(data);
    
    for jj = 1:n_series
        
        % waitbar
        if jj == 1
            f = waitbar((jj-1)/n_series, append('Processing image ', string(jj), ' of ', string(n_series), '...'));
        else
            waitbar((jj-1)/n_series, f, append('Processing image ', string(jj), ' of ', string(n_series), '...'));
        end
        
        series = data{jj,1};
        series_names = series(:,2);
        
        % find the images corresponding to each channel
        if n_channels == 3
            idx_405 = find(~cellfun(@isempty,(strfind(series_names, 'C=1/3'))));
            idx_488 = find(~cellfun(@isempty,(strfind(series_names, 'C=2/3'))));
            idx_555 = find(~cellfun(@isempty,(strfind(series_names, 'C=3/3'))));
        elseif n_channels == 4
            idx_405 = find(~cellfun(@isempty,(strfind(series_names, 'C=1/4'))));
            idx_488 = find(~cellfun(@isempty,(strfind(series_names, 'C=2/4'))));
            idx_555 = find(~cellfun(@isempty,(strfind(series_names, 'C=3/4'))));
            idx_647 = find(~cellfun(@isempty,(strfind(series_names, 'C=4/4'))));
        end
        if isempty(idx_647)
            error('Selected number of channels not detected. Check to make sure n_channels is correct.')
        end
        
        % move the images for each channel into their own stack
        im_stack_405 = series(idx_405,1);
        im_stack_488 = series(idx_488,1);
        im_stack_555 = series(idx_555,1);
        
        % convert from cell to 3D matrix
        im_stack_405 = cat(3,im_stack_405{:});
        im_stack_488 = cat(3,im_stack_488{:});
        im_stack_555 = cat(3,im_stack_555{:});
        
        % make MIP
        mip_405 = max(im_stack_405, [], 3);
        mip_488 = max(im_stack_488, [], 3);
        mip_555 = max(im_stack_555, [], 3);
        
        if n_channels == 4
            im_stack_647 = series(idx_647,1);
            im_stack_647 = cat(3,im_stack_647{:});
            mip_647 = max(im_stack_647, [], 3);
        end
        
        % store images in MATLAB
        index = length(ims_405) + 1;
        ims_405{index} = mip_405;
        ims_488{index} = mip_488;
        ims_555{index} = mip_555;
        if n_channels == 4
            ims_647{index} = mip_647;
        end
        
        % write MIP to folder
        % determine name of experiment (i.e name of folder)
        idx = regexp(path_to_file{1}, '/');
        experiment_name = extractAfter(path_to_file{1}, idx(end));
        experiment_name = extractBefore(experiment_name, filetype);
        
        folder_path = extractBefore(path_to_file{1}, experiment_name);
        cd(folder_path)
        
        if ~isfolder(experiment_name)
            mkdir(experiment_name)
        end
        
        cd(experiment_name)
        im_name = append('max_', experiment_name, '_n', string(jj));
        
        % write MIPs to folder
        imwrite(mip_405, append(im_name, '_405.tiff'))
        imwrite(mip_488, append(im_name, '_488.tiff'))
        imwrite(mip_555, append(im_name, '_555.tiff'))
        if n_channels == 4
            imwrite(mip_647, append(im_name, '_647.tiff'));
        end
    end
    waitbar(1, f, 'Done!');
    pause(3)
    close(f)
end

