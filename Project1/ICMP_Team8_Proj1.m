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
figure
set(gcf, 'color', [0.9412 0.9412 0.9412])
set(gca, 'color', [0.9412 0.9412 0.9412])
subplot(2,1,1)
plot(sample_pts, abp_hr10)
hold on 
% Onsets 
onsets = abp_hr10(onset_times);
plot(onset_times, onsets, '*', 'MarkerSize',10,'LineWidth', 1.5)

%  end of systole estimated wtih beat periods
EOF_bts_times = abpfeature_hr10(:,9);
EOF_bts = abp_hr10(EOF_bts_times);
plot(EOF_bts_times, EOF_bts, 'X', 'MarkerSize',10,'LineWidth', 1.5)

% end of systole estimated wtih lowest nonnegative slope
EOF_lns_times = abpfeature_hr10(:,11);
end_lns= abp_hr10(EOF_lns_times);
plot(EOF_lns_times, end_lns, 'O' , 'MarkerSize',10,'LineWidth', 1.5)

xlabel('Sample Number (125 samples = 1 s)');
ylabel('Arterial Blood Pressure, (mmHg)');
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
plot(onset_times_hr11, onsets_11, '*', 'MarkerSize',10,'LineWidth', 1.5)


EOF_bts_times_hr11 = abpfeature_hr11(:,9);
end_bts_11 = abp_hr11(EOF_bts_times_hr11);
plot(EOF_bts_times_hr11, end_bts_11, 'X', 'MarkerSize',10,'LineWidth', 1.5)


EOS_lns_times_hr11= abpfeature_hr11(:,11);
end_lns_hr11 = abp_hr11(EOS_lns_times_hr11);
plot(EOS_lns_times_hr11, end_lns_hr11, 'O' , 'MarkerSize',10,'LineWidth', 1.5)

xlabel('Sample Number (125 samples = 1 s)');
ylabel('Arterial Blood Pressure (mmHg)');
title('20 Pulses Starting at Hour 11')
% legend("ABP Waveform", "Onset", "End of Systole (0.3 * sqrt(Beat Period))", "End of Systole (lowest nonnegative slope)")

axis([0 1920 40 180])

hold off

%% Question 2
%% For patient 427
% Load in files
% ABP_File_427 = readtable('s00427-2977-04-30-13-51_ABP.txt');
% numeric_tbl_427 = readtable('s00427-2977-04-30-13-51n.txt');
% 
% % Grabbed indexes of points starting from hr 9 that would roughly give about 20 pulses (125 samples = 1 s)
% startIdx_427 = find(ABP_File_427.Var1 >= 32400, 1, 'first'); % 9 hours = 9*3600 s
% endIdx_427 = startIdx_26 + 1980;
% 
% ABP_subset_427 = ABP_File_427(startIdx_427:endIdx_427, :);
% abp_427 = ABP_subset_427.Var2;
% 
% % Obtained onset times using wabp function
% onset_times_427 = wabp(abp_427);
% 
% % Obtained ABP waveform features using abpfeature function
% abpfeature_427 = abpfeature(abp_427, onset_times_427);
% 
% % Number of sample pts for ~20 pulses
% sample_pts_427 = linspace(0, 1980, 1980);
% 
% figure
% subplot(2,1,1)
% plot(sample_pts_427, abp_427)
% hold on
% % Onsets
% onsets_427 = abp_427(onset_times_427);
% plot(onset_times_427, onsets_427, '*', 'MarkerSize', 10)
% 
% % End of systole estimated with beat periods
% EOF_bts_times_427 = abpfeature_427(:,9);
% EOF_bts_427 = abp_427(EOF_bts_times_427);
% plot(EOF_bts_times_427, EOF_bts_427, 'X', 'MarkerSize', 10)
% 
% % End of systole estimated with lowest nonnegative slope
% EOF_lns_times_427 = abpfeature_427(:,11);
% end_lns_427 = abp_427(EOF_lns_times_427);
% plot(EOF_lns_times_427, end_lns_427, 'O', 'MarkerSize', 10)
% 
% xlabel('Sample Number (125 samples = 1 s)')
% ylabel('Arterial Blood Pressure (mmHg)')
% title('Patient 427: 20 Pulses Starting at Hour 9')
% legend("ABP Waveform", "Onset", "End of Systole (0.3 * sqrt(Beat Period))", "End of Systole (lowest nonnegative slope)")
% axis([0 1980 40 180])
% hold off
% 
% 
% 
% %% For patient 20
% % Load in Files for patient 20
% ABP_File_20 = readtable('s00020-2567-03-30-17-47_ABP.txt');
% numeric_tbl_20 = readtable('s00020-2567-03-30-17-47n.txt');
% 
% % Choose a different hour window to show waveform variation
% startIdx_20 = find(ABP_File_20.Var1 >= 25200, 1, 'first');
% endIdx_20 = startIdx_20 + 1980;
% 
% ABP_subset_20 = ABP_File_20(startIdx_20:endIdx_20, :);
% abp_20 = ABP_subset_20.Var2;
% 
% % Onset times and features
% onset_times_20 = wabp(abp_20);
% abpfeature_20 = abpfeature(abp_20, onset_times_20);
% 
% sample_pts_20 = linspace(0, 1980, 1980);
% 
% subplot(2,1,2)
% plot(sample_pts_20, abp_20)
% hold on
% onsets_20 = abp_20(onset_times_20);
% plot(onset_times_20, onsets_20, '*', 'MarkerSize', 10)
% 
% EOF_bts_times_20 = abpfeature_20(:,9);
% EOF_bts_20 = abp_20(EOF_bts_times_20);
% plot(EOF_bts_times_20, EOF_bts_20, 'X', 'MarkerSize', 10)
% 
% EOF_lns_times_20 = abpfeature_20(:,11);
% end_lns_20 = abp_20(EOF_lns_times_20);
% plot(EOF_lns_times_20, end_lns_20, 'O', 'MarkerSize', 10)
% 
% xlabel('Sample Number (125 samples = 1 s)')
% ylabel('Arterial Blood Pressure (mmHg)')
% title('Patient 20: 20 Pulses Starting at Hour 7')
% legend("ABP Waveform", "Onset", "End of Systole (0.3 * sqrt(Beat Period))", "End of Systole (lowest nonnegative slope)")
% axis([0 1980 40 180])
% hold off



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
[co_5, to, told, fea] = estimateCO_v3(t_on,feat,beatq,estID,filt_order);


% calculated calibration C2 using first valid reference CO from 
% patient 20's numeric file (which was at 0.2 hrs)
C2 = numeric_tbl.CO(13)./co_5(1);

%Used C2 to calibrate estimated CO
CO_calibrated_5 = C2.*co_5;

% converted time from minutes to hours
time_hr_5 = to./60;

CO = numeric_tbl.CO;

% find  CO ouputs between hours 0 and 12 hours in numeric file
idx = find(CO ~= 0);
CO_time = numeric_tbl.ElapsedTime(idx) / 3600;
CO_val = CO(idx);

keep = CO_time <= 12;
CO_time = CO_time(keep);
CO_val = CO_val(keep);

% Performance of Algorithm
Pred_CO_5= interp1(time_hr_5, CO_calibrated_5, CO_time, 'nearest', 'extrap');
RMSE_5 = sqrt(mean((CO_val - Pred_CO_5).^2));
MAE_5 =  mean(abs(CO_val - Pred_CO_5));
r_5 = corr(CO_val, Pred_CO_5);

diffs_5 = CO_val - Pred_CO_5;
% Bias (mean difference)
bias_5 = mean(diffs_5);
% Standard deviation of differences
sd_diff_5= std(diffs_5);
% 95% Limits of Agreement
LoA_lower_5 = bias_5 - 1.96 * sd_diff_5;
LoA_upper_5 = bias_5 + 1.96 * sd_diff_5;



% Directional Agreement
delta_act_5 = diff(CO_val);       % reference changes
delta_pred_5  = diff(Pred_CO_5);   % predicted changes

dir_act_5 = sign(delta_act_5);
dir_pred_5   = sign(delta_pred_5);

% ignore zero changes in reference
nonzero_idx_5 = dir_act_5 ~= 0;
matches_5 = dir_act_5(nonzero_idx_5) == dir_pred_5(nonzero_idx_5);
dir_agr_5 = mean(matches_5);




PP = fea(:,5); % Pulse Pressure
PP_Num = interp1(time_hr_5, PP, CO_time, 'Linear', 'extrap');
PP_time = CO_time;


MAP = fea(:,6); % Mean Arterial Pressure
MAP_Num = interp1(time_hr_5, MAP, CO_time, 'Linear', 'extrap');
MAP_time = CO_time;


HR = 60*125./fea(:,7); % Heart Rate
HR_Num = interp1(time_hr_5, HR, CO_time, 'Linear', 'extrap');
HR_time = CO_time;


% Reproduction of Figure 4
figure
set(gcf, 'color', [0.9412 0.9412 0.9412])
set(gca, 'color', [0.9412 0.9412 0.9412])
% estimated CO plot
subplot(4,1,1)
plot(time_hr_5, CO_calibrated_5)
hold on
stem(CO_time, CO_val, 'MarkerFaceColor',[1, 0.5, 0], 'LineWidth', 2, 'Color',[1, 0.5, 0])
ylabel('CO (L/min)')
xlim([0 12])
ylim([0,10])
legend('Calibrated CO', "True CO")
title("Estimate of Continuous CO, PP, MAP, HR Using Liljestrand Estimator")
hold off 

% Pulse Pressure (PP) plot
subplot(4,1,2)
hold on
plot(time_hr_5, PP)
stem(PP_time, PP_Num, 'MarkerFaceColor',[1, 0.5, 0], 'LineWidth', 2, 'Color',[1, 0.5, 0])
ylabel('PP (mmHg)')
xlim([0 12])
ylim([0,120])
hold off

% mean arterial pressure (MAP) plot
subplot(4,1,3)
hold on
plot(time_hr_5, MAP)
stem(MAP_time, MAP_Num, 'MarkerFaceColor',[1, 0.5, 0], 'LineWidth', 2, 'Color',[1, 0.5, 0])
ylabel('MAP (mmHg)')
xlim([0 12])
hold off


% heart rate (HR) plot
subplot(4,1,4)
hold on 
plot(time_hr_5, HR)
stem(HR_time, HR_Num, 'MarkerFaceColor',[1, 0.5, 0], 'LineWidth', 2, 'Color',[1, 0.5, 0])
xlabel('Time [hours]');
ylabel('HR (bpm)');
xlim([0 12])
hold off 

%% Question 4

% Reproduction of Figure 4
figure
% estimated CO plot
subplot(4,1,1)
plot(time_hr_5, CO_calibrated_5)
hold on
stem(CO_time, CO_val, 'MarkerFaceColor',[1, 0.5, 0])
ylabel('CO (L/min)')
xlim([0 12])
ylim([2,10])
legend('Calibrated CO', "True CO")
title("Estimate of Continuous CO, PP, MAP, HR Using Liljestrand Estimator")
hold off 


% Estimator 2
estID = 2;
filt_order = 1;

% obtained estimated CO using Windkessel 1st order estimator
[co_2, to_2, told, fea] = estimateCO_v3(t_on,feat,beatq,estID,filt_order);


% calculated calibration C2 using first valid reference CO from 
% patient 20's numeric file (which was at 0.2 hrs)
C2 = numeric_tbl.CO(13)./co_2(1);

%Used C2 to calibrate estimated CO
CO_calibrated_2 = C2.*co_2;

% converted time from minutes to hours
time_hr_2 = to_2./60;
CO = numeric_tbl.CO;

% find  CO ouputs between hours 0 and 12 hours in numeric file
idx = find(CO ~= 0);
CO_time = numeric_tbl.ElapsedTime(idx) / 3600;
CO_val = CO(idx);

keep = CO_time <= 12;
CO_time = CO_time(keep);
CO_val = CO_val(keep);

% Performance of Algorithm
Pred_CO_2= interp1(time_hr_2, CO_calibrated_2, CO_time, 'nearest', 'extrap');
RMSE_2 = sqrt(mean((CO_val - Pred_CO_2).^2));
MAE_2 =  mean(abs(CO_val - Pred_CO_2));
r_2 = corr(CO_val, Pred_CO_2);


diffs_2 = CO_val - Pred_CO_2;
% Bias (mean difference)
bias_2 = mean(diffs_2);
% Standard deviation of differences
sd_diff_2= std(diffs_2);
% 95% Limits of Agreement
LoA_lower_2 = bias_2 - 1.96 * sd_diff_2;
LoA_upper_2 = bias_2 + 1.96 * sd_diff_2;


% Directional Agreement
delta_act_2 = diff(CO_val);       % reference changes
delta_pred_2  = diff(Pred_CO_2);   % predicted changes

dir_act_2 = sign(delta_act_2);
dir_pred_2   = sign(delta_pred_2);

% ignore zero changes in reference
nonzero_idx_2 = dir_act_2 ~= 0;
matches_2 = dir_act_2(nonzero_idx_2) == dir_pred_2(nonzero_idx_2);
dir_agr_2 = mean(matches_2);




% Reproduction of Figure 4
subplot(4,1,2)
% estimated CO plot
plot(time_hr_2, CO_calibrated_2)
hold on
stem(CO_time, CO_val, 'MarkerFaceColor',[1, 0.5, 0])
ylabel('CO (L/min)')
xlim([0 12])
ylim([2,10])
legend('Calibrated CO', "True CO")
title("Estimate of Continuous CO Using  Windkessel Estimator")
hold off 


estID = 7;
filt_order = 1;
% obtained estimated CO using Wesseling estimator
[co_7, to_7, told, fea] = estimateCO_v3(t_on,feat,beatq,estID,filt_order);


% calculated calibration C2 using first valid reference CO from 
% patient 20's numeric file (which was at 0.2 hrs)
C2 = numeric_tbl.CO(13)./co_7(1);

%Used C2 to calibrate estimated CO
CO_calibrated_7 = C2.*co_7;

% converted time from minutes to hours
time_hr_7 = to_7./60;

CO = numeric_tbl.CO;

% find CO ouputs between hours 0 and 12 hours in numeric file
idx = find(CO ~= 0);
CO_time = numeric_tbl.ElapsedTime(idx) / 3600;
CO_val = CO(idx);

keep = CO_time <= 12;
CO_time = CO_time(keep);
CO_val = CO_val(keep);

% Performance of Algorithm
Pred_CO_7= interp1(time_hr_7, CO_calibrated_7, CO_time, 'nearest', 'extrap');
RMSE_7 = sqrt(mean((CO_val - Pred_CO_7).^2));
MAE_7 =  mean(abs(CO_val - Pred_CO_7));
r_7 = corr(CO_val, Pred_CO_7);

diffs_7 = CO_val - Pred_CO_7;
% Bias (mean difference)
bias_7 = mean(diffs_7);
% Standard deviation of differences
sd_diff_7= std(diffs_7);
% 95% Limits of Agreement
LoA_lower_7 = bias_7 - 1.96 * sd_diff_7;
LoA_upper_7 = bias_7 + 1.96 * sd_diff_7;

% Directional Agreement
delta_act_7 = diff(CO_val);       % reference changes
delta_pred_7  = diff(Pred_CO_7);   % predicted changes

dir_act_7 = sign(delta_act_7);
dir_pred_7   = sign(delta_pred_7);

nonzero_idx_7 = dir_act_7 ~= 0;
matches_7 = dir_act_7(nonzero_idx_7) == dir_pred_7(nonzero_idx_7);
dir_agr_7 = mean(matches_7);




% Reproduction of Figure 4
subplot(4,1,3)
% estimated CO plot
plot(time_hr_7, CO_calibrated_7)
hold on
stem(CO_time, CO_val, 'MarkerFaceColor',[1, 0.5, 0])
ylabel('CO (L/min)')
xlim([0 12])
ylim([2,10])
legend('Calibrated CO', "True CO")
title("Estimate of Continuous CO Using Wesseling Estimator")
hold off 


subplot(4,1,4)
% estimated CO plot with all the different estimators
plot(time_hr_5, CO_calibrated_5)
hold on
plot(time_hr_7, CO_calibrated_7)
plot(time_hr_2, CO_calibrated_2)

stem(CO_time, CO_val, 'MarkerFaceColor','r')
ylabel('CO (L/min)')
xlim([0 12])
ylim([2,10])
legend('Calibrated CO', "True CO")
title("Estimate of Continuous CO, PP, MAP, HR Using Liljestrand Estimator")
hold off

figure 
set(gcf, 'color', [0.9412 0.9412 0.9412])
set(gca, 'color', [0.9412 0.9412 0.9412])
plot(time_hr_7, CO_calibrated_7, 'Color', [0.13, 0.39, 0.11]) % green
hold on
plot(time_hr_5, CO_calibrated_5, 'Color',[0.07,0.44,0.75]) % blue

plot(time_hr_2, CO_calibrated_2, 'Color', [1.00,0.91,0.39]) % yellow


stem(CO_time, CO_val, 'MarkerFaceColor',[1, 0.5, 0], 'MarkerEdgeColor', [1, 0.5, 0], 'LineWidth', 2, 'Color', [1, 0.5, 0])
ylabel('CO (L/min)')
xlabel('time [hours]')
xlim([0 12])
ylim([0,10])
legend('Wesseling (IC)', 'Liljestrand','Windkessel 1st order', "True CO")
% legend('Calibrated CO', "True CO")
title("Estimate of Continuous CO, PP, MAP, HR Using Various Estimators")
hold off

%% Question 5 

% Estimator 14 --> parlikar
estID = 14;
filt_order = 1;

% obtained estimated CO using parlikar estimator
[co_14, to_14, told, fea] = estimateCO_v3(t_on,feat,beatq,estID,filt_order);

C2 = numeric_tbl.CO(13)./co_14(1);
%Used C2 to calibrate estimated CO
CO_calibrated_14 = C2.*co_14;

% converted time from minutes to hours
time_hr_14 = to_14./60;

CO = numeric_tbl.CO;

% find  CO ouputs between hours 0 and 12 hours in numeric file
idx = find(CO ~= 0);
CO_time = numeric_tbl.ElapsedTime(idx) / 3600;
CO_val = CO(idx);

keep = CO_time <= 12;
CO_time = CO_time(keep);
CO_val = CO_val(keep);

% Performance of Algorithm
Pred_CO_14= interp1(time_hr_14, CO_calibrated_14, CO_time, 'nearest', 'extrap');
RMSE_14 = sqrt(mean((CO_val - Pred_CO_14).^2));
MAE_14 =  mean(abs(CO_val - Pred_CO_14));
r_14 = corr(CO_val, Pred_CO_14);

diffs_14 = CO_val - Pred_CO_14;
% Bias (mean difference)
bias_14 = mean(diffs_14);
% Standard deviation of differences
sd_diff_14= std(diffs_14);
% 95% Limits of Agreement
LoA_lower_14 = bias_14 - 1.96 * sd_diff_14;
LoA_upper_14 = bias_14 + 1.96 * sd_diff_14;


% Directional Agreement 
delta_act_14 = diff(CO_val);       % reference changes
delta_pred_14  = diff(Pred_CO_14);   % predicted changes

dir_act_14 = sign(delta_act_14);
dir_pred_14   = sign(delta_pred_14);

% ignore zero changes in reference
nonzero_idx_14 = dir_act_14 ~= 0;
matches_14 = dir_act_14(nonzero_idx_14) == dir_pred_14(nonzero_idx_14);
dir_agr_14 = mean(matches_14);




figure
% estimated CO plot Liljestrand C2 calibration
plot(time_hr_14, CO_calibrated_14, 'Color',[0.82, 0.02, 0.55])
hold on
plot(time_hr_5, CO_calibrated_5,'Color', [0.07,0.44,0.75])

stem(CO_time, CO_val, 'MarkerFaceColor',[1, 0.5, 0], 'MarkerEdgeColor', [1, 0.5, 0], 'LineWidth', 2, 'Color', [1, 0.5, 0])
ylabel('CO (L/min)')
xlim([0 12])
ylim([2,10])
legend('Parlikar', "Lijestrand", "True CO")
title("Estimate of Continuous CO, Parlikar vs Liljestrand")
hold off 



figure
% estimated CO plot parlikar C2 calibration
set(gcf, 'color', [0.9412 0.9412 0.9412])
set(gca, 'color', [0.9412 0.9412 0.9412])
plot(time_hr_2, CO_calibrated_2, 'Color', [1.00,0.91,0.39]) % yellow
hold on
plot(time_hr_7, CO_calibrated_7, 'Color', [0.13, 0.39, 0.11]) % green
plot(time_hr_5, CO_calibrated_5,'Color', [0.07,0.44,0.75])

plot(time_hr_14, CO_calibrated_14, 'Color',[0.82, 0.02, 0.55])

stem(CO_time, CO_val, 'MarkerFaceColor',[1, 0.5, 0], 'MarkerEdgeColor', [1, 0.5, 0], 'LineWidth', 2, 'Color', [1, 0.5, 0])
ylabel('CO (L/min)')
xlabel('time [hours]')

xlim([0 12])
ylim([2,10])
legend('Windkessel 1st order','Wesseling (IC)', "Lijestrand",'Parlikar', "True CO")
title("Estimate of Continuous CO, Parlikar vs Other Algorithms")
hold off 


co_14_x = interp1(time_hr_14, co_14, CO_time, 'nearest', 'extrap');
co_5_x = interp1(time_hr_5, co_5, CO_time, 'nearest', 'extrap');


% Calibration testing
C1_14 = (CO_val.' * co_14_x) ./ (co_14_x.' * co_14_x);
CO_C1_calibrated_14 = C1_14.*co_14;

C1_5 = (CO_val.' * co_5_x) ./ (co_5_x.' * co_5_x);

CO_C1_calibrated_5 = C1_5.*co_5;



% patient 20
% Performance of Algorithm C1 Parlikar
Pred_CO_14_C1= interp1(time_hr_14, CO_C1_calibrated_14, CO_time, 'nearest', 'extrap');
RMSE_14_C1 = sqrt(mean((CO_val - Pred_CO_14_C1).^2));
MAE_14_C1 =  mean(abs(CO_val - Pred_CO_14_C1));
r_14_C1 = corr(CO_val, Pred_CO_14_C1);

% diffs_14_C1 = CO_val - Pred_CO_14_C1;
% % Bias (mean difference)
% bias_14_C1 = mean(diffs_14_C1);
% % Standard deviation of differences
% sd_diff_14_C1= std(diffs_14_C1);
% % 95% Limits of Agreement
% LoA_lower_14_C1 = bias_14_C1 - 1.96 * sd_diff_14_C1
% LoA_upper_14_C1 = bias_14_C1 + 1.96 * sd_diff_14_C1
% 

% delta_act_14_C1 = diff(CO_val);       % reference changes
% delta_pred_14_C1  = diff(Pred_CO_14);   % predicted changes
% 
% dir_act_14_C1 = sign(delta_act_14_C1);
% dir_pred_14_C1   = sign(delta_pred_14_C1);
% 
% % ignore zero changes in reference
% nonzero_idx_14_C1 = dir_act_14_C1 ~= 0;
% matches_14_C1 = dir_act_14_C1(nonzero_idx_14_C1) == dir_pred_14_C1(nonzero_idx_14_C1);
% dir_agr_14_C1 = mean(matches_14_C1)


% Performance of Algorithm C1 Liljestrand
Pred_CO_5_C1= interp1(time_hr_5, CO_C1_calibrated_5, CO_time, 'nearest', 'extrap');
RMSE_5_C1 = sqrt(mean((CO_val - Pred_CO_5_C1).^2));
MAE_5_C1 =  mean(abs(CO_val - Pred_CO_5_C1));
r_5_C1 = corr(CO_val, Pred_CO_5_C1);

% diffs_5_C1 = CO_val - Pred_CO_5_C1;
% % Bias (mean difference)
% bias_5_C1 = mean(diffs_5_C1);
% % Standard deviation of differences
% sd_diff_5_C1= std(diffs_14_C1);
% % 95% Limits of Agreement
% LoA_lower_5_C1 = bias_5_C1 - 1.96 * sd_diff_5_C1
% LoA_upper_5_C1 = bias_5_C1 + 1.96 * sd_diff_5_C1
% 
% 
% 
% delta_act_5_C1 = diff(CO_val);       % reference changes
% delta_pred_5_C1  = diff(Pred_CO_14);   % predicted changes
% 
% dir_act_5_C1 = sign(delta_act_5_C1);
% dir_pred_5_C1   = sign(delta_pred_5_C1);
% 
% % ignore zero changes in reference
% nonzero_idx_5_C1 = dir_act_5_C1 ~= 0;
% matches_5_C1 = dir_act_5_C1(nonzero_idx_5_C1) == dir_pred_5_C1(nonzero_idx_5_C1);
% dir_agr_5_C1 = mean(matches_5_C1)


% Calibration comparisions between Liljestrand C1 vs C2
figure
set(gcf, 'color', [0.9412 0.9412 0.9412])
set(gca, 'color', [0.9412 0.9412 0.9412])
subplot(2,1,1)
plot(time_hr_5, CO_calibrated_5,'Color',[0.07,0.44,0.75])
hold on
plot(time_hr_5, CO_C1_calibrated_5,'Color', [1.00, 0, 0])

stem(CO_time, CO_val, 'MarkerFaceColor',[1, 0.5, 0], 'MarkerEdgeColor', [1, 0.5, 0], 'LineWidth', 2, 'Color', [1, 0.5, 0])
ylabel('CO (L/min)')
xlabel('time [hours]')
xlim([0 12])
ylim([0,10])
legend('C2','C1',"True CO")
title("Estimate of Continuous CO Liljestrand, C2 vs C1")
hold off 

% Calibration comparisions between Parlikar C1 vs C2

figure
set(gcf, 'color', [0.9412 0.9412 0.9412])
set(gca, 'color', [0.9412 0.9412 0.9412])
subplot(2,1,1)
plot(time_hr_14, CO_calibrated_14,'Color',[0.82, 0.02, 0.55])
hold on
plot(time_hr_14, CO_C1_calibrated_14,'Color', [0.15, 0.15, 0.15])

stem(CO_time, CO_val, 'MarkerFaceColor',[1, 0.5, 0], 'MarkerEdgeColor', [1, 0.5, 0], 'LineWidth', 2, 'Color', [1, 0.5, 0])
ylabel('CO (L/min)')
xlabel('time [hours]')
xlim([0 12])
ylim([0,10])
legend('C2','C1',"True CO")
title("Estimate of Continuous CO Parlikar, C2 vs C1")
hold off 



myFlag=1;   % Program finished


return
