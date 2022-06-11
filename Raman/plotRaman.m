%Plot Raman Spectrum

lam_laser = 660E-9; % wavelength of excitation laser in m

dirname = "20220524_Raman/";
[x0, I0]=  readvars(dirname + "20220524_paracetamol_300mW.csv");
[x1, I1] = readvars(dirname + "20220524-Aline_test_BsubPill.csv");
[x2, I2] = readvars(dirname + "20220524-Aline_test_BsubPill2.csv");
[x3, I3] = readvars(dirname + "20220524-Aline_test_BsubPill3.csv");
[x4, I4] = readvars(dirname + "20220524-Aline_test_BsubPill4.csv");
[x5, I5] = readvars(dirname + "20220524-Aline_test_BsubPill5.csv");

figure; hold on
I_adj_1 = (I1-min(I1))/(max(I1)-min(I1));
I_adj_2 = (I2-min(I2))/(max(I2)-min(I2));
I_adj_3 = (I3-min(I3))/(max(I3)-min(I3));
I_adj_4 = (I4-min(I4))/(max(I4)-min(I4));
I_adj_5 = (I5-min(I5))/(max(I5)-min(I5));


I_pill_mean = (I_adj_1+I_adj_2+I_adj_3+I_adj_4+I_adj_5)/5;
plot(x1, I_adj_1, x2, I_adj_2, x3, I_adj_3, x4, I_adj_4, x5, I_adj_5);

plot(I_pill_mean);
legend('1','2','3','4','5', 'mean');





figure; 
hold on


[Base1, Corrected_Spectrum1]=baseline(I1);
plot(x1,I1,'b-',x1,Base1,'r--',x1,Corrected_Spectrum1,'g-.');

[Base2, Corrected_Spectrum2]=baseline(I2);
plot(x2,I2,'b-',x1,Base2,'r--',x1,Corrected_Spectrum2,'g-.');

[Base3, Corrected_Spectrum3]=baseline(I3);
plot(x3,I3,'b-',x1,Base3,'r--',x1,Corrected_Spectrum3,'g-.');

[Base4, Corrected_Spectrum4]=baseline(I4);
plot(x4,I4,'b-',x1,Base4,'r--',x1,Corrected_Spectrum4,'g-.');

[Base5, Corrected_Spectrum5]=baseline(I5);
plot(x5,I5,'b-',x1,Base5,'r--',x1,Corrected_Spectrum5,'g-.');

figure;hold on;
[Base, Corrected_Spectrum]=baseline(I0);
plot(x0,I0,'b-',x0,Base,'r--',x0,Corrected_Spectrum,'g');



%% Display the output of BEADS
fc = 0.015;     % fc : cut-off frequency (cycles/sample)
d = 1;          % d : filter order parameter (d = 1 or 2)
r = 8;          % r : asymmetry parameter

% Regularization parameters
amp = 0.1;
lam0 = 0.5*amp;
lam1 = 5*amp;
lam2 = 4*amp;
% I0(1:25) = []; % chop regions around zero which might be an unaccurate offset
I0 = I0-min(I0);
N = length(I0);

wavelengths = (lam_laser)+(0:1/7000:(N-1)/7000)*10e-9;
wavenumbers = 1/lam_laser - 1./wavelengths;

[X0, f0, cost] = beads(I0-min(I0),d,fc, r, lam0, lam1, lam2);

ylim1 = [-20 inf+20];
xlim1 = [0 N+20];

figure(1); hold on
clf

subplot(4, 1, 1)
plot(I0)
title('Data');
xlim(xlim1)
ylim(ylim1)

subplot(4, 1, 2)
plot(I0,'color', [1 1 1]*0.7)
line(1:N, f0, 'LineWidth', 1)
legend('Data', 'Baseline')
legend boxoff
title(['Baseline, as estimated by BEADS', ' (r = ', num2str(r), ', fc = ', num2str(fc), ', d = ', num2str(d),')'])
xlim(xlim1)
ylim(ylim1)

subplot(4, 1, 3)
plot(X0)
title('Baseline-corrected data')
xlim(xlim1)
ylim(ylim1)
% set(gca,'ytick', ylim1)


subplot(4, 1, 4)
plot(I0 - X0 - f0)
title('Residual')
xlim(xlim1)
ylim(ylim1)
% set(gca,'ytick', ylim1)

orient tall
print -dpdf example

%% Compare results

N = length(I0);
fc = 0.006;     % fc : cut-off frequency (cycles/sample)
d = 1;          % d : filter order parameter (d = 1 or 2)
r = 6;          % r : asymmetry parameter
% Regularization parameters
amp = 0.8;
lam0 = 0.5*amp;
lam1 = 5*amp;
lam2 = 4*amp;

wavelengths = (lam_laser)+(0:1/7000:(N-1)/7000)*10e-9;
wavenumbers = 1/lam_laser - 1./wavelengths;

[X0, f0, cost] = beads(I0-min(I0),d,fc, r, lam0, lam1, lam2);

close figure 4; figure(4);
h4 = axes;
set(h4, 'Xdir', 'reverse');

[Base0, Corrected_Spectrum0]=baseline(I0);

plot(wavenumbers, Corrected_Spectrum0); hold on
plot(wavenumbers, X0, 'LineWidth', 1);
legend('Corrected Spectrm (baseline fn)', 'Corrected value (BEADS)');
% legend boxoff
title(['Baseline Corrected Data, BEADS parameters: ', ' (r = ', num2str(r), ', fc = ', num2str(fc), ', d = ', num2str(d),')'])
xlabel("Wavenumber (cm^{-1})");
ylabel("Intensity");


%%
figure; hold on;
[Base, Corrected_Spectrum]=baseline(I_pill_mean);
plot(x1,I_pill_mean,'b-',x1,Base,'r--',x1,Corrected_Spectrum,'g');
%% Display the output of BEADS
fc = 0.025;     % fc : cut-off frequency (cycles/sample)
d = 1;          % d : filter order parameter (d = 1 or 2)
r = 10;          % r : asymmetry parameter

% Regularization parameters
amp = 0.1;
lam0 = 0.5*amp;
lam1 = 5*amp;
lam2 = 4*amp;
% I1(1:30) = []; % chop regions around zero which might be an unaccurate offset
% I1 = I1-min(I1);
N = length(I1);
[X1, f1, cost] = beads(I1-min(I1),d,fc, r, lam0, lam1, lam2);

ylim1 = [-50 inf];
xlim1 = [0 N+50];

figure(6); hold on
clf

subplot(4, 1, 1)
plot(I1-min(I1))
title('Data');
xlim(xlim1)
ylim(ylim1)

subplot(4, 1, 2)
plot(I1-min(I1),'color', [1 1 1]*0.7)
line(1:N, f1, 'LineWidth', 1)
legend('Data', 'Baseline')
legend boxoff
title(['Baseline, as estimated by BEADS', ' (r = ', num2str(r), ', fc = ', num2str(fc), ', d = ', num2str(d),')'])
xlim(xlim1)
ylim(ylim1)

subplot(4, 1, 3)
plot(X1)
title('Baseline-corrected data')
xlim(xlim1)
ylim(ylim1)
set(gca,'ytick', ylim1)


subplot(4, 1, 4)
plot(I1-min(I1) - X1 - f1)
title('Residual')
xlim(xlim1)
ylim(ylim1)
set(gca,'ytick', ylim1)

orient tall
print -dpdf example
