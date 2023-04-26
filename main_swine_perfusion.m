clear all
close all
clc                   

%% ---------- Global settings
dsa_dirpath = 'dsa';
dir_mask_artery_before_mesh = 'mask_aif';
dir_masks = 'masks';
result_dir = 'results\exp_opus87';
frame_rate = 15;
if ~exist(result_dir, 'dir')
   mkdir(result_dir)
end
intermediate_result_mat_dir = fullfile(result_dir, 'mat');
if ~exist(intermediate_result_mat_dir, 'dir')
   mkdir(intermediate_result_mat_dir)
end
acquisition_info_csv = 'acquisitions.csv';
tt = readtable(acquisition_info_csv);  
dsa_filelist = dir(fullfile(dsa_dirpath, '**\*.IMA'));
results_mat = 'results';
landmarks = {'artery', 'parenchymal', 'vein'};

%% ----------- Generating parametric images
disp("Generating parameter images ...")
for dd = 1:length(dsa_filelist)
    tic;
    dsa_filepath = fullfile(dsa_filelist(dd).folder, dsa_filelist(dd).name);
    dsa_seq = squeeze(dicomread(dsa_filepath));
    dsa_seq = dsa_seq(:, :, 6:end); % removing first five frames to avoid mask frame just in case; assuming the last dim is time dim.
    dsa_seq(dsa_seq==0) = median(dsa_seq(:));
    ds_info = dicominfo(dsa_filepath);
    
    dsa_seq = double(255 * mat2gray(dsa_seq));
    dsa_seq = 255 - dsa_seq;
    dsa_seq = permute(dsa_seq, [3 1 2]);  % move time axis to the 1st dim.
    
    acquisition_info = tt(strcmp(tt.Patient, ds_info.PatientName.FamilyName) & tt.SeriesNumber==ds_info.SeriesNumber, :);
    result = table2struct(acquisition_info);
    result.filestem = dsa_filelist(dd).name(1:end-4);
    result.MinIP = squeeze(min(255 - dsa_seq));
    result.SeqID = strcat(ds_info.PatientName.FamilyName, "-", num2str(ds_info.SeriesNumber, '%02d'));
    result.Patient = ds_info.PatientName.FamilyName;
    disp(result.SeqID)
    
    % AIF
    AIF_mask_filepath = fullfile(dir_aif_mask, strcat(dsa_filelist(dd).name(1:end-4), '.png'));
    result.aif_mask = imread(AIF_mask_filepath) > 0; 

    %subtraction of baseline signal
    dsa_seq=dsa_seq-repmat(median(dsa_seq(:)),[size(dsa_seq,1) 1 1]);
    dsa_seq(dsa_seq<0) = 0;

    % apply random injection delay into the DSA series 
    dsa_seq = introduce_injection_delay(dsa_seq, result.aif_mask, round(acquisition_info.InjectionDelay * frame_rate));

    [result.CBF, result.CBV, result.Tmax, result.MTT, tics, irfs] = modelfree_deconvolution(dsa_seq, result.aif_mask, 'dt', 1/frame_rate, 'prefilter', 1);

    % aif characteristics
    result.aif = mean(tics(:,result.aif_mask==1),2);
    result.aif_volume = sum(result.aif);
    [result.aif_peak, result.aif_ttp] = max(result.aif);
    result.aif_ttp = result.aif_ttp/frame_rate;

    % tic derived parameter maps
    result.volume = squeeze(sum(tics));
    [result.peak, result.ttp] = max(tics); 
    result.peak = squeeze(result.peak);
    result.ttp = squeeze(result.ttp/frame_rate);
    result.intensity = squeeze(mean(tics));

    save(fullfile(intermediate_result_mat_dir, strcat(result.SeqID, '_tics.mat')), 'tics');
    save(fullfile(intermediate_result_mat_dir, strcat(result.SeqID, '_irfs.mat')), 'irfs');
    
    toc 
end
% save(results_mat, 'results');
disp("---------------------------------------------")

%% Save the average TICs and IRFs of masked ROIs to results (time consuming at loading the mat files)
for rr = 1:numel(results)
    tic;
    for landmark = landmarks
        % load saved tics and irfs
        landmark = char(landmark);
        disp(strcat("-----------", results(rr).Patient, "-", num2str(results(rr).ID), "-", landmark))
        dsa_filepath = fullfile(dsa_filelist(rr).folder, dsa_filelist(rr).name);
        ds_info = dicominfo(dsa_filepath);
        landmark_mask_filepath = fullfile(dir_masks, landmark, strcat(results(rr).filestem, '.png'));
        results(rr).(strcat(landmark, '_mask')) = imread(landmark_mask_filepath) > 0;
        
        load(fullfile(result_dir, 'mat', strcat(ds_info.PatientName.FamilyName, "-", num2str(ds_info.SeriesNumber, '%02d'), '_tics.mat')));
        results(rr).(strcat(landmark, '_tic')) = mean(tics(:,results(rr).(strcat(landmark, '_mask'))==1),2);

        load(fullfile(result_dir, 'mat', strcat(ds_info.PatientName.FamilyName, "-", num2str(ds_info.SeriesNumber, '%02d'), '_irfs.mat')));
        results(rr).(strcat(landmark, '_irf')) = mean(irfs(:,results(rr).(strcat(landmark, '_mask'))==1),2);
    end
    toc;
    save(results_mat, 'results');
end

%% Post-processing of parameter maps
disp("Post-processing parameter maps...");
for rr = 1:numel(results)
    CBF_threshold = 0.02;
    results(rr).CBF(results(rr).CBV<=0) = 0;
    results(rr).Tmax(results(rr).CBV<=0) = 0;      
    results(rr).MTT(results(rr).CBV<=0) = 0;
    results(rr).CBV(results(rr).CBV<=0) = 0;

    results(rr).MTT(results(rr).CBF<=CBF_threshold) = 0;
    results(rr).Tmax(results(rr).CBF<=CBF_threshold) = 0;

    results(rr).MTT = medfilt2(results(rr).MTT);
    results(rr).Tmax = medfilt2(results(rr).Tmax);
end

%% Gather parameter min and max values
disp("Calculating min and max values of parameter maps...");
minmax = struct('CBV_max', -inf, 'Tmax_max', -inf, 'MTT_max', -inf, ...
    'CBV_min', inf, 'CBF_min', inf, 'Tmax_min', inf, 'MTT_min', inf);

for result = results
    % gather and update min max values
    minmax.CBV_max = max(minmax.CBV_max, max(result.CBV, [], 'all'));
    minmax.CBF_max = max(minmax.CBF_max, max(result.CBF, [], 'all'));
    minmax.Tmax_max = max(minmax.Tmax_max, max(result.Tmax, [], 'all'));
    minmax.MTT_max = max(minmax.MTT_max, max(result.MTT, [], 'all'));

    minmax.CBV_min = min(minmax.CBV_min, min(result.CBV, [], 'all'));
    minmax.CBF_min = min(minmax.CBF_min, min(result.CBF, [], 'all'));
    minmax.Tmax_min = min(minmax.Tmax_min, min(result.Tmax, [], 'all'));
    minmax.MTT_min = min(minmax.MTT_min, min(result.MTT, [], 'all')); 
end
disp("done");

%% Save parametric image visualizations
disp("Saving parametric image visualizations ...")
for result = results
    disp(strcat("-----------", result.Patient, "-", num2str(result.ID)))
    figure('Visible', 'off');
    t = tiledlayout(1, 5, 'Padding', 'none', 'TileSpacing', 'none'); xlabel(t, result.SeqID);
    MinIP = imoverlay(uint8(result.MinIP),result.aif_mask, [0.8500 0.3250 0.0980]); % orange
    nexttile; imshow(MinIP); title('MinIP'); xticks([]);yticks([]);
    nexttile; image(uint8(255 * mat2gray(result.CBV, [minmax.CBV_min minmax.CBV_max]))); title('CBV'); axis image; xticks([]); yticks([]);
    nexttile; image(uint8(255 * mat2gray(result.CBF, [minmax.CBF_min minmax.CBF_max]))); title('CBF'); axis image; xticks([]); yticks([]);
    nexttile; image(uint8(255 * mat2gray(result.Tmax, [minmax.Tmax_min minmax.Tmax_max]))); title('Tmax'); axis image; xticks([]); yticks([]);
    nexttile; image(uint8(255 * mat2gray(result.MTT, [minmax.MTT_min minmax.MTT_max]))); title('MTT'); axis image; xticks([]); yticks([]);
    exportgraphics(t, fullfile(result_dir, strcat(result.SeqID, '.png')), 'Resolution', 1200);
    % --------------------------------
end
disp("done");

%% Compute avg params over ROIs and save to CSV.
disp("Calculating average CBF, CBV, Tmax, MTT, volume, peak, and ttp, over masked ROIs ...")    
for rr = 1:numel(results)
    for landmark = landmarks
        landmark = char(landmark);
        disp(strcat("-----------", results(rr).Patient, "-", num2str(results(rr).ID), "-", landmark))

        results(rr).(strcat(landmark, '_volume')) = mean(results(rr).volume(results(rr).(strcat(landmark, '_mask'))==1));
        results(rr).(strcat(landmark, '_peak')) = mean(results(rr).peak(results(rr).(strcat(landmark, '_mask'))==1));
        results(rr).(strcat(landmark, '_ttp')) = mean(results(rr).ttp(results(rr).(strcat(landmark, '_mask'))==1));
        for param = {'CBV', 'CBF', 'Tmax', 'MTT'}
            results(rr).(strcat(landmark, '_', char(param))) = mean(results(rr).(char(param))(results(rr).(strcat(landmark, '_mask'))==1));  % average param values of masked ROI
        end
        results(rr).(strcat(landmark, '_intensity')) = mean(results(rr).intensity(results(rr).(strcat(landmark, '_mask'))==1));
    end
end
result_table = struct2table(results);
disp(result_table)
non_writeable_columns = {'MinIP', 'artery_mask', 'parenchymal_mask', 'vein_mask'};
non_writeable_columns = cat(2, non_writeable_columns, {'aif_mask', 'aif', ...
    'CBV', 'CBF', 'MTT', 'Tmax', 'volume', 'peak', 'ttp', 'intensity', ...
    'artery_tic', 'parenchymal_tic', 'vein_tic', 'artery_irf', 'parenchymal_irf', 'vein_irf'});
writetable(removevars(result_table, non_writeable_columns), fullfile(result_dir, 'results_deconv_avg.csv'))