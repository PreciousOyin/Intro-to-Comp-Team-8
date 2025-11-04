function myFlag = ICMP_Team8_Proj1()
clear all
close all
myFlag = 0;

% Load in Files
ABP_File = readtable('s00020-2567-03-30-17-47_ABP.txt');
numeric_tbl = readtable('s00020-2567-03-30-17-47n.txt');

%% Question 1

% Grabbed indexes of points starting from hr 10 that would roughly
% give about 20 pulses

% Var1 is time in seconds and Var2 is ABP
startIdx = find(ABP_File.Var1 == 36000, 1, 'first');
endIdx = startIdx + 2000;

ABP_subset = ABP_File(startIdx:endIdx, :);
abp_hr10 =ABP_subset.Var2;


% Obtained Onset Times using wabp function
onset_times = wabp(abp_hr10);


% Obtained ABP waveform features using abpfeature function
abpfeature_hr10 = abpfeature(abp_hr10, onset_times);


% Number of samples pts for the 20 pulses
sample_pts = linspace(0, 2001, 2001);

subplot(2,1,1)
plot(sample_pts, abp_hr10)
hold on 
% Onsets 
onsets = abp_hr10(onset_times);
plot(onset_times, onsets, '*', 'MarkerSize',10)

%  end of systole estimated wtih beat periods
EOF_bts_times = abpfeature_hr10(:,9);
EOF_bts = abp_hr10(EOF_bts_times);
plot(EOF_bts_times, EOF_bts, 'X', 'MarkerSize',10)

% end of systole estimated wtih lowest nonnegative slope
EOF_lns_times = abpfeature_hr10(:,11);
end_lns= abp_hr10(EOF_lns_times);
plot(EOF_lns_times, end_lns, 'O' , 'MarkerSize',10)

xlabel('Sample Number');
ylabel('Arterial Blood Pressure');
title('Hour 10 Pulses')

axis([0 2100 40 180]);

hold off


% Repeated same steps for hr 11

% hour 11
startIdx = find(ABP_File.Var1 == 39600, 1, 'first');
end_idx = startIdx + 2000;


ABP_subset_11 = ABP_File(startIdx:end_idx, :);

abp_hr11 =ABP_subset_11.Var2;


onset_times_hr11 = wabp(abp_hr11);

abpfeature_hr11 = abpfeature(abp_hr11, onset_times_hr11);


sample_pts = linspace(0, 2001, 2001);

subplot(2,1,2)
plot(sample_pts, abp_hr11)
hold on 
onsets_11 = abp_hr11(onset_times_hr11);
plot(onset_times_hr11, onsets_11, '*', 'MarkerSize',10)


EOF_bts_times_hr11 = abpfeature_hr11(:,9);
end_bts_11 = abp_hr11(EOF_bts_times_hr11);
plot(EOF_bts_times_hr11, end_bts_11, 'X', 'MarkerSize',10)


EOS_lns_times_hr11= abpfeature_hr11(:,11);
end_lns_hr11 = abp_hr11(EOS_lns_times_hr11);
plot(EOS_lns_times_hr11, end_lns_hr11, 'O' , 'MarkerSize',10)

xlabel('Sample Number');
ylabel('Arterial Blood Pressure');
title('Hour 11 Pulses')

axis([0 2100 40 180])
hold off





%% Question 3

%data of ABP for the first 12 hours
startIdx = 1;
end_idx = find(ABP_File.Var1 == 43200, 1, 'first');
abp_subset = ABP_File(startIdx:end_idx, :);

% Var1 is time in seconds and Var2 is ABP
abp = abp_subset.Var2;

% onsets times and abp features for 12 hr time frame
onset_12hr = wabp(abp);
abpfeature_12hr = abpfeature(abp, onset_12hr);
t_on = onset_12hr;
feat  = abpfeature_12hr;

% using jSQI function to get SQI of each beat (BeatQ_12hr) and other 
% features for TCO measurements (PP, MAP, HR)
[BeatQ_12hr, frac_good_bts]= jSQI(abpfeature_12hr, onset_12hr, abp);
beatq = BeatQ_12hr;

% Liljestrand PP/(Psys+Pdias) estimator
estID = 5;

filt_order = 1;

% obtained estimated CO using Liljestrand estimator
[co, to, told, fea] = estimateCO_v3(t_on,feat,beatq,estID,filt_order);


% calculated calibration C3 using first valid reference CO from 
% patient 20's numeric file (which was at hr 12)
C2 = numeric_tbl.CO(13)./co(1);

%Used C2 to calibrate estimated CO
CO_calibrated = C2.*co;

% converted time from minutes to hours
time_hr = to./60;


% size(CO_calibrated)
% size(time_hr)
% size(fea)


% Reproduction of Figure 4
figure
% estimated CO plot
subplot(4,1,1)
plot(time_hr, CO_calibrated)
ylabel('CO')
ylim([2,10])

% estimated Pulse Pressure (PP) plot
subplot(4,1,2)
plot(time_hr, fea(:,5))
ylabel('PP')
ylim([20,120])

% mean arterial pressure (MAP) plot
subplot(4,1,3)
plot(time_hr, fea(:,6))
ylabel('MAP')

% heart rate (HR) plot
subplot(4,1,4)
plot(time_hr, 60*125./fea(:,6))
xlabel('Time [hours]');
ylabel('HR');


myFlag=1;   % Program finished


return