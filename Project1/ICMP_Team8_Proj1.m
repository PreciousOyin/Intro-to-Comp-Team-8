% Team Name: Team 8
% Team Members: Precious Oyinloye, Anna Vargas, Armaan Bilimoria, Lei Yan
% CPA
% File Number 1 of 1

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
startIdx = find(ABP_File.Var1 >= 36000, 1, 'first');
endIdx = startIdx + 1979;

ABP_subset = ABP_File(startIdx:endIdx, :);
abp_hr10 =ABP_subset.Var2;


% Obtained Onset Times using wabp function
onset_times = wabp(abp_hr10);


% Obtained ABP waveform features using abpfeature function
abpfeature_hr10 = abpfeature(abp_hr10, onset_times);


% Number of samples pts for the 20 pulses
sample_pts = linspace(0, 1980, 1980);

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

xlabel('Sample Number (125 samples = 1 s)');
ylabel('Arterial Blood Pressure');
title('20 Pulses Starting at Hour 10')

axis([0 1940 40 180]);
legend("ABP Waveform", "Onset", "End of Systole (0.3 * sqrt(Beat Period))", "End of Systole (lowest nonnegative slope)")
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

xlabel('Sample Number (125 samples = 1 s)');
ylabel('Arterial Blood Pressure');
title('20 Pulses Starting at Hour 11')
legend("ABP Waveform", "Onset", "End of Systole (0.3 * sqrt(Beat Period))", "End of Systole (lowest nonnegative slope)")

axis([0 1920 40 180])
hold off

%% Question 2
%% For patient 427
% Load in files
ABP_File_427 = readtable('s00427-2977-04-30-13-51_ABP.txt');
numeric_tbl_427 = readtable('s00427-2977-04-30-13-51n.txt');

% Grabbed indexes of points starting from hr 9 that would roughly give about 20 pulses (125 samples = 1 s)
startIdx_427 = find(ABP_File_427.Var1 >= 32400, 1, 'first'); % 9 hours = 9*3600 s
endIdx_427 = startIdx_26 + 1980;

ABP_subset_427 = ABP_File_427(startIdx_427:endIdx_427, :);
abp_427 = ABP_subset_427.Var2;

% Obtained onset times using wabp function
onset_times_427 = wabp(abp_427);

% Obtained ABP waveform features using abpfeature function
abpfeature_427 = abpfeature(abp_427, onset_times_427);

% Number of sample pts for ~20 pulses
sample_pts_427 = linspace(0, 1980, 1980);

figure
subplot(2,1,1)
plot(sample_pts_427, abp_427)
hold on
% Onsets
onsets_427 = abp_427(onset_times_427);
plot(onset_times_427, onsets_427, '*', 'MarkerSize', 10)

% End of systole estimated with beat periods
EOF_bts_times_427 = abpfeature_427(:,9);
EOF_bts_427 = abp_427(EOF_bts_times_427);
plot(EOF_bts_times_427, EOF_bts_427, 'X', 'MarkerSize', 10)

% End of systole estimated with lowest nonnegative slope
EOF_lns_times_427 = abpfeature_427(:,11);
end_lns_427 = abp_427(EOF_lns_times_427);
plot(EOF_lns_times_427, end_lns_427, 'O', 'MarkerSize', 10)

xlabel('Sample Number (125 samples = 1 s)')
ylabel('Arterial Blood Pressure (mmHg)')
title('Patient 427: 20 Pulses Starting at Hour 9')
legend("ABP Waveform", "Onset", ...
    "End of Systole (0.3 * sqrt(Beat Period))", ...
    "End of Systole (lowest nonnegative slope)")
axis([0 1980 40 180])
hold off



%% For patient 20
% Load in Files for patient 20
ABP_File_20 = readtable('s00020-2567-03-30-17-47_ABP.txt');
numeric_tbl_20 = readtable('s00020-2567-03-30-17-47n.txt');

% Choose a different hour window to show waveform variation
startIdx_20 = find(ABP_File_20.Var1 >= 25200, 1, 'first');
endIdx_20 = startIdx_20 + 1980;

ABP_subset_20 = ABP_File_20(startIdx_20:endIdx_20, :);
abp_20 = ABP_subset_20.Var2;

% Onset times and features
onset_times_20 = wabp(abp_20);
abpfeature_20 = abpfeature(abp_20, onset_times_20);

sample_pts_20 = linspace(0, 1980, 1980);

subplot(2,1,2)
plot(sample_pts_20, abp_20)
hold on
onsets_20 = abp_20(onset_times_20);
plot(onset_times_20, onsets_20, '*', 'MarkerSize', 10)

EOF_bts_times_20 = abpfeature_20(:,9);
EOF_bts_20 = abp_20(EOF_bts_times_20);
plot(EOF_bts_times_20, EOF_bts_20, 'X', 'MarkerSize', 10)

EOF_lns_times_20 = abpfeature_20(:,11);
end_lns_20 = abp_20(EOF_lns_times_20);
plot(EOF_lns_times_20, end_lns_20, 'O', 'MarkerSize', 10)

xlabel('Sample Number (125 samples = 1 s)')
ylabel('Arterial Blood Pressure (mmHg)')
title('Patient 20: 20 Pulses Starting at Hour 7')
legend("ABP Waveform", "Onset", ...
    "End of Systole (0.3 * sqrt(Beat Period))", ...
    "End of Systole (lowest nonnegative slope)")
axis([0 1980 40 180])
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


% calculated calibration C2 using first valid reference CO from 
% patient 20's numeric file (which was at 0.2 hrs)
C2 = numeric_tbl.CO(13)./co(1);

%Used C2 to calibrate estimated CO
CO_calibrated = C2.*co;

% converted time from minutes to hours
time_hr = to./60;

CO = numeric_tbl.CO;

% find  CO ouputs between hours 0 and 12 hours in numeric file
idx = find(CO ~= 0);
CO_time = numeric_tbl.ElapsedTime(idx) / 3600;
CO_val = CO(idx);

keep = CO_time <= 12;
CO_time = CO_time(keep);
CO_val = CO_val(keep);


PP = fea(:,5); % Pulse Pressure
MAP = fea(:,6); % Mean Arterial Pressure
HR = 60*125./fea(:,7); % Heart Rate



% Reproduction of Figure 4
figure
% estimated CO plot
subplot(4,1,1)
plot(time_hr, CO_calibrated)
hold on
stem(CO_time, CO_val, 'MarkerFaceColor',[1, 0.5, 0])
ylabel('CO')
xlim([0 12])
ylim([2,10])
legend('Calibrated CO', "True CO")
title("Estimate of Continuous CO, PP, MAP, HR Using Liljestrand Estimator")
hold off 

% Pulse Pressure (PP) plot
subplot(4,1,2)
plot(time_hr, PP)
ylabel('PP (mmHg)')
ylim([20,120])

% mean arterial pressure (MAP) plot
subplot(4,1,3)
plot(time_hr, MAP)
ylabel('MAP (mmHg)')

% heart rate (HR) plot
subplot(4,1,4)
plot(time_hr, HR)
xlabel('Time [hours]');
ylabel('HR (bpm)');


myFlag=1;   % Program finished


return
