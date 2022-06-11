function obj = findpeaks_sort(data, N_max_peaks, fake_peaks_width)
% findpeaks Find local peaks in data, sort them and extract maximal peaks
% INPUT: 
    %   data: DATA (1D array)
    %   N_max_peaks: NUMBER OF PEAKS to select
    %   fake_peaks_width: WIDTH below which a peak is considered too thin
        %   to be an independent peak and not a bump/ noise/ resonnance
        %   onto a larger peak
% OUTPUT
    %   obj- struct with fields: (see findpeaks() function)
        %   pks: peak amplitudes
        %   locs: locations (indexes) of the peaks (ascending)
        %   w: width of the peaks
        %   p: prominence of the peaks
        %   maxpeaks: N_max_peaks long, only maximal
        %   maxlocs: corresponding locations (index)
        %   w_max: corresponding peak widths
        %   p_max: corresponding peak prominence

if ~exist('fake_peaks_width', 'var') || isempty(fake_peaks_width)
    fake_peaks_width = 1;
end

[pks, locs, w, p] = findpeaks(data);
            % can precise sampling rate if needed
            % but can also keep simple index and relate manually to wavelength

if ~exist('N_max_peaks', 'var') || isempty(N_max_peaks) % check if we precised how many peaks we want
    N_max_peaks = length(pks);
end


%store peaks (all, sorted in order of index)
obj.pks = pks;
obj.locs = locs; % index of peak --> transform to wavelength
obj.w = w;
obj.p = p;


% sort by highest intensity
[pks_sorted, I] = sort(pks, 'descend');
maxlocs = locs(I);
w_max  = w(I);
p_max = p(I);

maxpeaks = pks_sorted(1:N_max_peaks); % take only the first N peaks (defined at the start)
maxlocs = maxlocs(1:N_max_peaks);

[maxlocs, I] = sort(maxlocs, 'ascend');
maxpeaks = maxpeaks(I);
w_max  = w_max(I);
p_max = p_max(I);

% remove "bumps" fake peaks
I = find(w_max < fake_peaks_width);
maxpeaks(I) = [];
maxlocs(I) = [];
w_max(I) = [];
p_max(I) =[];

obj.maxpks = maxpeaks;
obj.maxlocs = maxlocs; % index of peak --> transform to wavelength
obj.w_max = w_max; % width of the peaks
obj.p_max = p_max;


end