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





%% Question 3 (IDK rn)

startIdx = 1;
end_idx = find(ABP_File.Var1 == 43200, 1, 'first');

T_subset = ABP_File(startIdx:end_idx, :);

abp_test =T_subset.Var2;

onset_test = wabp(abp_test);

abpfeature_test = abpfeature(abp_test, onset_test);


t_on = onset_test;
feat  = abpfeature_test;
[BeatQ_test, frac_good_bts]= jSQI(abpfeature_test, onset_test, abp_test);

beatq = BeatQ_test;

estID = 5;

filt_order = 1;

[co, to, told, fea] = estimateCO_v3(t_on,feat,beatq,estID,filt_order);

% numeric_tbl.CO(13)

C3 = numeric_tbl.CO(13)./co;

% CO_calibrated = C3.*co;



myFlag=1;   % Program finished


return