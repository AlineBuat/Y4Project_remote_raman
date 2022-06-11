%% Analysing Raman Spectra
%% Aline BUAT
%% 01/06/2022

clear;
close all;


% Requires the BEADS toolbox
% files should be .csv named starting with the date "YYYYMMDD_"
% file names should also contain the keywords to the label for later
% convenience

%addpath("BEADS_toolbox\");


%% load spectra from files
% select files
dirname = "20220601_Raman\";
file_extension = ".csv";
% filenames = dirname + "*" + file_extension;
filenames = path_name+ dirname + "*" + file_extension;
files = dir(filenames); % open files from directory
varnames = ["Pixel", "Intensity"];
files(end+1) = dir(path_name+"20220524_Raman\20220524_paracetamol_100mW.csv");
files(end+1) = dir(path_name+"20220524_Raman\20220524_paracetamol_300mW.csv");


% tiling of image
N_files=length(files);
if ~isprime(N_files)
    k_num = factor(N_files);
end

m_rows=floor(sqrt(N_files)*1.5); n_cols=ceil(sqrt(N_files)/1.5);
if m_rows*n_cols<N_files
    m_rows=m_rows+1; fprintf('Adjusted figure size hehe\n');
end

table_stacks = [];
intensity_stack = {};


% import table in table format from all files
for n_it = 1:N_files
    file = files(n_it); % file is a local copy
    filename = file.name;
    % label string by type of sample
    if contains(filename, "b_sub", 'IgnoreCase', true) || contains(filename, "bsub", 'IgnoreCase',true)
        files(n_it).label = "b_sub";
        if contains(filename, "optical_glue", 'IgnoreCase', true) || contains(filename, "opticalglue", 'IgnoreCase', true)
            files(n_it).label2 = "optical_glue";
        end
    elseif contains(filename, "paracetamol", 'IgnoreCase',true)
        files(n_it).label = "paracetamol";
    elseif (contains(filename, "optical_glue", 'IgnoreCase', true) || contains(filename, "opticalglue", 'IgnoreCase', true))
        files(n_it).label = "optical_glue";
    else files(n_it).label = "unknown";
    end


    filefolder = file.folder;
    T = readtable(filefolder+"\"+filename);
    allVars = 1:width(T);
    T = renamevars(T, allVars, varnames);
    table_stacks(n_it).name = filename;
    table_stacks(n_it).T = T;
    files(n_it).spectra.intensity = T.Intensity;
    intensity_stack{n_it} = T.Intensity;

end

%% Parameters for baseline extraction
fc = 0.011;     % fc : cut-off frequency (cycles/sample)
d = 1;          % d : filter order parameter (d = 1 or 2)
% Positivity bias (peaks are positive)
r = 15;          % r : asymmetry parameter

% Regularization parameters
amp = 0.2;
lam0 = 0.5*amp;
lam1 = 5*amp;
lam2 = 4*amp;

% how many maximal peaks to take
N_max_peaks = 14;

N_max_peaks_2 = 20; % for general characterization

%initial delay to remove
N_delay = 27;
N_end = 500;

% paracetamol peaks (in m-1)
para_true_peaks = [1651, 1612, 1559]; % Raman shift, wavenumbers in cm-1
para_true_peaks = para_true_peaks(end:-1:1);
para_true_peaks_m = para_true_peaks*1E2; % in m-1
lam_laser = 660; %  laser wavelength in nm
lam_laser_m = lam_laser*1E-9; % laser wavelength in m

para_true_peaks_wl_m = 1./(1./lam_laser_m - para_true_peaks_m);
para_true_peaks_wl_nm = para_true_peaks_wl_m*1E9;


%% start having fun with the files
close(findobj('type','figure','name', 'Example_Initial_Data'))
figure('name','Example_Initial_Data');
plot(table_stacks(10).T.Intensity);

close(findobj('type','figure','name', 'Raw_Raman_traces'))
figure('name', 'Raw_Raman_traces');
for n_it=1:N_files
    subplot(m_rows, n_cols, n_it);
    plot(intensity_stack{n_it});
%     axis off;
    title(files(n_it).name(1:end-4),'Interpreter','none');
end

%% clean up the spectra

I_adj_stack = cell(size(intensity_stack));
cleaned_spectra =  cell(size(intensity_stack));
baselines =  cell(size(intensity_stack));
costs =  cell(size(intensity_stack));
cleaned_spectra2 =  cell(size(intensity_stack));
baselines2 =  cell(size(intensity_stack));

for n_it = 1:length(intensity_stack)
    

    I = intensity_stack{n_it};
    I(1:N_delay)=[];
    I_adj = (I-min(I));
    I_adj_stack{n_it} = I_adj;

    [X, f, cost] = beads(I_adj, d, fc, r, lam0, lam1, lam2);
    [Base, Corrected_Spectrum] = baseline(I_adj);
    cleaned_spectra{n_it} = X;
    baselines{n_it} = f;
    costs{n_it} = cost;
    cleaned_spectra2{n_it} = Corrected_Spectrum;
    baselines2{n_it} = Base;

    files(n_it).spectra.I_adj = I_adj;
    files(n_it).spectra.cleaned = X;
    files(n_it).spectra.baseline = f;
end


%% Plot cleaned up traces


% clean up process
close(findobj('type','figure','name', 'Raman_baseline'))
f_intensity_adj = figure('name', 'Raman_baseline'); hold off;
f_intensity_adj.Position = [15, 20, 1200, 750];
for n_it=1:length(intensity_stack)
    subplot(m_rows, n_cols, n_it); hold off
    plot(I_adj_stack{n_it}, 'LineWidth', 0.75); hold on;
    plot(baselines{n_it},'-.' , 'LineWidth', 0.75);
    plot(cleaned_spectra{n_it}, 'LineWidth', 0.75, 'Color', '#FF9933');
%     axis off;
    title(files(n_it).name(10:end-4),'Interpreter','none');
    xlabel('Pixel');
    ylabel('Intensity');
end
sgtitle("Baseline Substraction for Raman Spectra");

close(findobj('type','figure','name', 'Raman_cleaned'))
f_multiplot_cleaned = figure('name', 'Raman_cleaned');
% f_multiplot_cleaned.Position = [15, 20, 1200, 750];

% cleaned spectra only
for n_it=1:length(intensity_stack)
    subplot(m_rows, n_cols, n_it); hold off;
    plot(cleaned_spectra{n_it}); hold on;
%     axis off;
    title(files(n_it).name(10:end-4),'Interpreter','none');
end
sgtitle('Raman Spectra of samples for varying laser power and exposure time') 

%% plot with the second method
% f_intensity_adj2 = figure(5); hold off;
% for n_it=1:length(intensity_stack)
%     subplot(m, n, n_it); hold off
%     plot(I_adj_stack{n_it}); hold on;
%     plot(baselines2{n_it});
%     plot(cleaned_spectra2{n_it})
% %     axis off;
%     title(files(n_it).name(10:end-4),'Interpreter','none');
% end

% f_multiplot_cleaned2 = figure(6);
% f_multiplot_cleaned2.Position = [20, 20, 1200, 750];
% 
% for n_it=1:length(intensity_stack)
%     subplot(m, n, n_it); hold off;
%     plot(cleaned_spectra2{n_it}); hold on;
% %     axis off;
%     title(files(n_it).name(10:end-4),'Interpreter','none');
% end

%% CONCLUSION:
% the baseline function sometimes works, but not consistently on data.
% The beads function is somewhat efficient with the current parameters
% However higher harmonic frequencies appear in the spectrum for an unknown
% reason (maybe amplified by opticla glue + coverslip) and are picked up as Raman peaks




% %% find frequencies for the spectrum.
% % Hopefully the regular peaks should be associated with a high intensity
% % Fourier Transfer frequency band
% Fs = 1; % sampling frequency
% N = 1; % filter order
% N2 = 2;
% w0 = 0.047;
% bw =w0/10;
% band1 = [0.015, 0.03]; % frequencies to cut
% [num1, den1] = butter(N,band1/(Fs/2),'stop'); % Generate filter 
% % [num1, den1]= iirnotch(w0, bw);
% 
% Freqs = cell(size(I_adj_stack));
% for n_it = 1:N_files
% %     I = cleaned_spectra{n_it};
%     I = I_adj_stack{n_it};
%     L = length(I);
%     f_ax = [-pi+pi/L:2*pi/L:pi-pi/L];
%     %      n = 2^nextpow2(L);
%     F = fftshift(fft(I));
%     Freqs{n_it}.F = F;
%     Freqs{n_it}.f_ax = f_ax;
%     deriv_stack{n_it} = deriv;
%     I_f = filter(num1, den1, I);
%     F2 = fftshift(fft(I_f));
%     I_f_stack{n_it} = I_f;
% end
% 
% 
% hFreqs = figure(7); 
% subplot(2,2,1)
% % hFreqs.Position = [20, 20, 1200, 750];
% hold off;
% plot(cleaned_spectra{7}, 'LineWidth', 1); hold on;
% legend('Orignal')
% subplot(2,2,3)
% plot(abs(I_f_stack{7}), 'LineWidth', 1);
% legend('Filtered')
% subplot(1,2,2); hold off;
% plot(Freqs{7}.f_ax, log(abs(Freqs{7}.F)) ); hold on;
% plot(Freqs{7}.f_ax, log(abs(fftshift(fft(I_f_stack{7})))) );
% legend('Original', 'Filtered');
% 
% figure(8);
% for n_it=1:length(intensity_stack)
%     subplot(m_rows, n_cols, n_it); hold off;
%     plot(Freqs{n_it}.f_ax, log(abs(Freqs{n_it}.F))); hold on;
%     title(files(n_it).name(10:end-4),'Interpreter','none');
% end

% conclusion: while the frequencies can be cut, this heavily interferes
% with the peaks (suspected "TRUE") and shall not be continued.
% peak amplitude is suspected to sometimes resonate with the periodic
% noise, which lead to it being amplified

%% Finding peaks
Peaks(N_files).pks = [];
Peaks(N_files).locs = [];
Peaks(N_files).w =[];
Peaks(N_files).p=[];
Peaks(N_files).maxpeaks = [];
Peaks(N_files).maxlocs = [];


for n_it = 1:N_files
    data = cleaned_spectra{n_it}(1:end-N_end);
    [pks,locs,w,p] = findpeaks(data);
    [pks_sorted, I] = sort(pks, 'descend');
    maxlocs = locs(I);

    maxpeaks = pks_sorted(1:N_max_peaks); % take only the first N peaks (defined at the start)
    maxlocs = maxlocs(1:N_max_peaks);

    [maxlocs, I] = sort(maxlocs, 'ascend');
    maxpeaks = maxpeaks(I);
    Peaks(n_it).pks = pks;
    Peaks(n_it).locs = locs;
    Peaks(n_it).w = w;
    Peaks(n_it).p = p;

    Peaks(n_it).maxpeaks = maxpeaks; 
    Peaks(n_it).maxlocs = maxlocs;
    files(n_it).spectra.Peaks = Peaks(n_it);
end

%% plot maximal peak positions

close(findobj('type','figure','name', 'Peak_localizations'))
f_peaks_locs = figure('name', 'Peak_localizations');
for n_it=1:length(intensity_stack)
    subplot(m_rows, n_cols, n_it); hold off;
    plot(cleaned_spectra{n_it}(1:end-N_end)); hold on;
    ys = min(cleaned_spectra{n_it}):max(cleaned_spectra{n_it}(1:end-N_end))*1.1;

    for k_it = 1:length(Peaks(n_it).maxlocs)
        plot(Peaks(n_it).maxlocs(k_it)*ones(size(ys)), ys, '--r');
    end
%     axis off;
    title(files(n_it).name(10:end-4),'Interpreter','none');
    xlabel("Pixel");
    ylabel("Intensity");
end

sgtitle('Raman Spectra of samples for varying laser power and exposure time') 

%% Calibrate wavelength to wavenumber using paracetamol peaks

para_clean_spectrum = files(end).spectra.cleaned(1:end-N_end);
para_peaks = findpeaks_sort(para_clean_spectrum, 15, 2.5);



medium_peaks = para_peaks.maxlocs(find(para_peaks.maxlocs<800 & para_peaks.maxlocs>400)); % identify peak regions: 3 peaks of interest

as = (para_true_peaks_wl_nm' - lam_laser)./medium_peaks; % nm/pixel
alpha = mean(as);
% a = (para_true_peaks_wl_nm(end)-lam_laser)/468; % nm/pixel

pixels = 1:length(para_clean_spectrum);
wavelengths =  lam_laser + alpha*pixels;
wavenums_m = 1/(lam_laser*1E-9) - 1./(wavelengths*1E-9); % in m-1
wavenums_cm = wavenums_m*1E-2; % in cm-1

% find maximal peaks (adjusted)
max_peaks_wavenum = wavenums_cm(para_peaks.maxlocs);


% plots
close(findobj('type','figure','name', 'paracetamol_adjust_wavenums'));
figure('name', 'paracetamol_adjust_wavenums'); hold on
subplot(2,1,1); hold on;
plot(pixels, para_clean_spectrum);
ys = min(para_clean_spectrum):max(para_clean_spectrum)*1.1;

for k_it = 1:length(para_peaks.maxlocs)
        plot(para_peaks.maxlocs(k_it)*ones(size(ys)), ys, '--r');
end

ylabel("Intensity");
xlabel("Pixel");
xlim([375, 550]);
title("Peaks from Paracetamol Spectrum Used for calibration")

% subplot(2,2,2);hold on;
% plot(wavelengths, para_clean_spectrum);
% ylabel("Intensity");
% xlabel("approximate lambdas");

subplot(2,1,2); hold on;
plot(wavenums_cm, para_clean_spectrum);
ylabel("Intensity");
xlabel("Adjusted wavenumbers (cm^{-1})");
title("Adjusted Paracetamol Raman spectrum")


%% Attempt classification

% plot all peaks side to side

close(findobj('type','figure','name', 'peaks'))
figure('name', 'peaks');
hold on; grid on
for n_it = 1:N_files
    N_peak_possible = length(Peaks(n_it).pks);
%     ys = n_it*ones(1, N_max_peaks);
    ys = Peaks(n_it).pks;
    plot(Peaks(n_it).locs, ys, '.', 'LineWidth',1, 'DisplayName', files(n_it).label);
end
xlabel("Pixel");
ylabel("Peak intensity");


%% Identify paracetamol peaks and adjust for wavelength to pixel scaling
para_num = find([files.label] =="paracetamol");

para_files = files(para_num);
mean_para_spectrum =0;
close(findobj('type','figure','name', 'paracetamol_spectrum'))
f_fig = figure('name', 'paracetamol_spectrum'); hold on;
peak_colour = ['r', 'b', 'g', 'y', 'k'];
for n_it = 1:length(para_files)
    spectrum = para_files(n_it).spectra.cleaned(1:end-N_end);
    spectrum = spectrum/max(spectrum)*10;
    plot(spectrum, 'LineWidth', 1); hold on;
    
    mean_para_spectrum = mean_para_spectrum + spectrum;

    ys = min(spectrum):max(spectrum)*1.1;
    pks_locs = para_files(n_it).spectra.Peaks.maxlocs;
    for k_it = 1:length(Peaks(n_it).maxlocs)
        plot(pks_locs*ones(size(ys)), ys, '--', 'Color', peak_colour(n_it));
    end
end

mean_para_spectrum = mean_para_spectrum/max(mean_para_spectrum)*10;

%% Average all b_subtillis plots

b_sub_num = find([files.label] =="b_sub");

b_sub_files = files(b_sub_num);

close(findobj('type','figure','name', 'b_sub_spectrum_stacked'))
f_fig = figure('name', 'b_sub_spectrum_stacked'); hold on;
peak_colour = ['r', 'b', 'g'];


mean_b_sub_spectrum = 0;

for n_it = 1:length(b_sub_files)
    spectrum = b_sub_files(n_it).spectra.cleaned(1:end-N_end);
    spectrum = spectrum/max(spectrum(150:800))*10;
%     spectrum = spectrum/max(spectrum)*10;
    plot(spectrum, 'LineWidth', 1); hold on;
    mean_b_sub_spectrum = mean_b_sub_spectrum + spectrum;
    ys = min(spectrum):max(spectrum)*1.1;
    pks_locs = b_sub_files(n_it).spectra.Peaks.maxlocs;
    for k_it = 1:length(Peaks(n_it).maxlocs)
%         plot(pks_locs*ones(size(ys)), ys, '--', 'Color', peak_colour(n_it));
    end
    

end

% mean_b_sub_spectrum = mean_b_sub_spectrum/length(b_sub_files);
mean_b_sub_spectrum = mean_b_sub_spectrum/max(mean_b_sub_spectrum)*10;

%% locate the peaks (b subtillis pill)
b_sub = findpeaks_sort(mean_b_sub_spectrum, [], 3);
maxlocs = b_sub.maxlocs;

close(findobj('type','figure','name', 'mean_b_sub_spectrum'))
f_fig = figure('name', 'mean_b_sub_spectrum'); hold on;
plot(mean_b_sub_spectrum);


ys = min(mean_b_sub_spectrum):max(mean_b_sub_spectrum)*1.1;
for k_it = 1:length(maxlocs)
        plot(maxlocs*ones(size(ys)), ys, '--r');
end


%% average optical glue
opt_glue_num = find([files.label] =="optical_glue");

opt_glue_files = files(opt_glue_num);

close(findobj('type','figure','name', 'optical_glue_spectrum'))
f_fig = figure('name', 'optical_glue_spectrum'); hold on;
peak_colour = ['r', 'b', 'g'];


mean_opt_glue_spectrum = 0;

for n_it = 1:length(opt_glue_files)
    spectrum = opt_glue_files(n_it).spectra.cleaned(1:end-N_end);
    spectrum = spectrum/max(spectrum(150:800))*10;
%     spectrum = spectrum/max(spectrum)*10;
    
    mean_opt_glue_spectrum = mean_opt_glue_spectrum + spectrum;
    
    pixels = 1:length(spectrum);
    wavelengths =  lam_laser + alpha*pixels;
    wavenums_m = 1/(lam_laser*1E-9) - 1./(wavelengths*1E-9); % in m-1
    wavenums_cm = wavenums_m*1E-2; % in cm-1
    
    plot(wavenums_cm, spectrum, 'LineWidth', 1); hold on;
end
title("Optical Glue Spectrum");
xlabel("Wavenumbers (cm^{-1})");
ylabel("Adjusted Intensity")
legend;
mean_opt_glue_spectrum = mean_opt_glue_spectrum/length(opt_glue_files);
meanval_opt_glue = mean(mean_opt_glue_spectrum);

%% locate the peaks (optical glue)
opt_glue = findpeaks_sort(mean_opt_glue_spectrum, [], 3);
maxlocs = opt_glue.maxlocs;

close(findobj('type','figure','name', 'mean_opt_glue_spectrum'))
f_fig = figure('name', 'mean_opt_glue_spectrum'); hold on;
plot(mean_opt_glue_spectrum);


ys = min(mean_opt_glue_spectrum):max(mean_opt_glue_spectrum)*1.1;
for k_it = 1:length(maxlocs)
        plot(maxlocs*ones(size(ys)), ys, '--r');
end


%% Attempt to substract the background
x = [1:length(mean_opt_glue_spectrum)]';
y = x.*maxpeaks(end)/2/length(mean_opt_glue_spectrum);
% y = (x)*2.78/1056;
adj_mean_opt = mean_opt_glue_spectrum-y;
adj_end_max = maxpeaks(end) - y(maxlocs(end)); % = adj_mean_opt(maxlocs(end))
Xend = maxlocs(end);
X1 = 650;

a = log(adj_end_max)/(maxlocs(end)-X1);
% x2 = [X1:length(mean_opt_glue_spectrum)]';

my_exp = exp(a*(x-X1));
reverso_exp = my_exp(end:-1:1);
my_exp2 = exp(a/1.5*(x-X1));
reverso_exp2 = my_exp2(end:-1:1);

close(findobj('type','figure','name', 'adj_glue_spectrum'))
figure('name', 'adj_glue_spectrum'); 
hold on;
plot(mean_opt_glue_spectrum);
plot(y);
plot(adj_mean_opt);
plot(x, exp(a*(x-X1)), x, -exp(a*(x-X1)) );
plot(x, reverso_exp);
plot(x, adj_mean_opt.*reverso_exp, 'LineWidth', 1);


adj_glue = adj_mean_opt.*reverso_exp;
% adj_glue = adj_mean_opt.*reverso_exp+ adj_mean_opt/1.5;
adj_glue = adj_mean_opt.*reverso_exp2;
plot(x, adj_glue,'--', 'LineWidth', 1, 'Color', 'r');
adj_glue = adj_glue/max(adj_glue)*10; % adjust the range

legend;

%% 
% mean_cleaned_b_sub = mean_b_sub_spectrum - (mean_opt_glue_spectrum)/1.3;
mean_cleaned_b_sub = mean_b_sub_spectrum - adj_glue/1.5;
close(findobj('type','figure','name', 'corrected_spectrum_b_sub'))
figure('name', 'corrected_spectrum_b_sub');
plot(mean_cleaned_b_sub);

%% compare peak locations and extract non-glue peaks
maxlocs_b = b_sub.maxlocs;
maxlocs_o = opt_glue.maxlocs;
comparison_mat = -1000*ones(length(maxlocs_b), length(maxlocs_o));


for n_it = 1:length(maxlocs_b)
    for k_it = 1:length(maxlocs_o)
        comparison_mat(n_it, k_it) = maxlocs_b(n_it)-maxlocs_o(k_it);
    end
end

[row, col] = find(abs(comparison_mat) < 8);

% compare to make sure
non_unique_peaks_b = maxlocs_b(unique(row));
unique_peaks_b = maxlocs_b(setdiff(1:end, row));
non_unique_peaks_o = maxlocs_o(unique(col));
unique_peaks_o = maxlocs_o(setdiff(1:end, col));


%% plot unique peaks of b_sub
pixels = 1:length(mean_b_sub_spectrum);
wavelengths =  lam_laser + alpha*pixels;
wavenums_m = 1/(lam_laser*1E-9) - 1./(wavelengths*1E-9); % in m-1
wavenums_cm = wavenums_m*1E-2; % in cm-1


close(findobj('type','figure','name', 'b_sub_unique_peaks'))
f_fig = figure('name', 'b_sub_unique_peaks'); hold on;
plot(wavenums_cm, mean_b_sub_spectrum);


ys = min(mean_b_sub_spectrum):max(mean_b_sub_spectrum)*1.1;
for k_it = 1:length(unique_peaks_b)
        plot(wavenums_cm(unique_peaks_b(k_it))*ones(size(ys)), ys, '--r');
end
xlabel('Wavenumber (cm^{-1})');
ylabel('Intensity (adjusted)');
title('Plot of Mean Sample (pill+glue) Signal', "Red: Unique peaks distinct from glue+glass slide");

%% unique peaks for all acquired spectra
locs_o = opt_glue.locs;

for i_it = 1:length(b_sub_files)
    Peaks = b_sub_files(i_it).spectra.Peaks;
    locs_b = Peaks.locs; % look for all peaks :(
    
    for n_it = 1:length(locs_b)
        for k_it = 1:length(locs_o)
            comparison_mat(n_it, k_it) = locs_b(n_it)-locs_o(k_it);
        end
    end

    u_peaks_b{i_it} = locs_b(setdiff(1:end, row));
    u_peaks_o{i_it} = locs_o(setdiff(1:end, col));

end

close(findobj('type','figure','name', 'b_sub_unique_peaks_all'))
f_fig = figure('name', 'b_sub_unique_peaks_all'); hold on;

N_b_sub_files = length(b_sub_files);
m_b_sub_rows =floor(sqrt(N_b_sub_files)*1.5); n_b_sub_cols=ceil(sqrt(N_b_sub_files)/1.5);
if m_b_sub_rows*n_b_sub_cols<N_b_sub_files
    m_rows=m_rows+1;
end



for n_it=1:length(b_sub_files)
    subplot(m_b_sub_rows, n_b_sub_cols, n_it); hold off
    spectrum = b_sub_files(n_it).spectra.cleaned;
    plot(spectrum, 'LineWidth', 0.75); hold on;

    ys = min(spectrum):max(spectrum)*1.1;
    for k_it = 1:length(u_peaks_b{n_it})
        plot(u_peaks_b{n_it}(k_it)*ones(size(ys)), ys, '--r');
    end
    
    title(b_sub_files(n_it).name(10:end-4),'Interpreter','none');
    xlabel('Wavenumber (cm^{-1})');
    ylabel('Intensity (adjusted)');
end
sgtitle("Unique Peaks of B. sub Raman Spectra");



%% plot with wavenumbers instead of wavelengths


close(findobj('type','figure','name', 'Raman_baseline'))
f_intensity_adj = figure('name', 'Raman_baseline'); hold off;
f_intensity_adj.Position = [15, 20, 1200, 750];
for n_it=1:length(intensity_stack)
    subplot(m_rows, n_cols, n_it); hold off
    pixels = 1:length(I_adj_stack{n_it});
    wavelengths =  lam_laser + alpha*pixels;
    wavenums_m = 1/(lam_laser*1E-9) - 1./(wavelengths*1E-9); % in m-1
    wavenums_cm = wavenums_m*1E-2; % in cm-1
    plot(wavenums_cm, I_adj_stack{n_it}, 'LineWidth', 0.75); hold on;
    plot(wavenums_cm, baselines{n_it},'-.' , 'LineWidth', 0.75);
    plot(wavenums_cm, cleaned_spectra{n_it}, 'LineWidth', 0.75, 'Color', '#FF9933');
%     axis off;
    title(files(n_it).name(10:end-4),'Interpreter','none');
    xlabel('Wavenumber (cm^{-1})');
    ylabel('Intensity');
end
sgtitle("Baseline Substraction for Raman Spectra");


%% Baseline extraction Example for paracetamol

close(findobj('type','figure','name', 'Raman_baseline_example'))
f_intensity_adj = figure('name', 'Raman_baseline_example'); hold on;
I_adj = para_files(end).spectra.I_adj;
baseline = para_files(end).spectra.baseline;
cleaned = para_files(end).spectra.cleaned;
peaks_locs = para_files(end).spectra.Peaks.maxlocs;

pixels = 1:length(I_adj);
wavelengths =  lam_laser + alpha*pixels;
wavenums_m = 1/(lam_laser*1E-9) - 1./(wavelengths*1E-9); % in m-1
wavenums_cm = wavenums_m*1E-2; % in cm-1

plot(wavenums_cm, I_adj, 'LineWidth', 1);
plot(wavenums_cm, baseline, 'LineWidth', 1, "Color", '#7E2F8E');
plot(wavenums_cm, cleaned, 'LineWidth', 1, "Color", '#D95319');

ys = min(I_adj)-max(I_adj)*0.025:max(I_adj)*1.1;
for k_it = 1:length(peaks_locs)
    plot(wavenums_cm(peaks_locs(k_it))*ones(size(ys)), ys, ':', "Color", '#A2142F');
end

title("Cleanup and peaks of Raman spectrum of Paracetamol")

xlabel("Wavenumbers (cm^{-1})");
ylabel("Intensity");
legend("Original", "Baseline", "Corrected", 'Peaks')