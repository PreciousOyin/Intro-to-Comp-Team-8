function myFlag = ICMP_Team8_Proj3()

clear all;
close all;

% Driver: Simplified model for Stefanini et al
% one isoform of VEGF; one receptor; one antibody
% three tissues (blood (central), tumor, rest of body)
% 2025-11-18

% index array (for readability of equations)
n.VEGF_b = 1;
n.VEGFR2_b = 2;
n.VEGFVEGFR2_b = 3;
n.Ab_b = 4;
n.VEGFAb_b = 5;
n.VEGF_t = 6;
n.VEGFR2_t = 7;
n.VEGFVEGFR2_t = 8;
n.Ab_t = 9;
n.VEGFAb_t = 10;
n.VEGF_r = 11;
n.VEGFR2_r = 12;
n.VEGFVEGFR2_r = 13;
n.Ab_r = 14;
n.VEGFAb_r = 15;

% parameters - compartment volumes (specifically, interstitial space)
p.Vol_b = 5*.6; % 5 litres, 60% available
p.Vol_t = 1*.611; % 1 litre, 61% available
p.Vol_r = 40*.0816; % 40 litres, 8% available

% parameters - clearance and internalization
p.kcl_V  = 0.0648; % min-1
p.kcl_VA = 2.2e-5;
p.kcl_A  = 2.2e-5;

p.kintR2  = 0.0168; %min-1
p.kintVR2 = 0.0168;

% parameters - production
p.VEGFprod_b = 0;
p.VEGFprod_t = 3; 
p.VEGFprod_r = 10; 

p.VEGFR2prod_b = 0; 
p.VEGFR2prod_t = 2.85e2*p.kintR2; % p.VEGFR2prod_t = goal conc * p.kintR2
p.VEGFR2prod_r = 2.2e3*p.kintR2; % p.VEGFR2prod_n = goal conc * p.kintR2

% parameters - binding
p.konVR  = 6e-4; % pM-1 min-1 
p.koffVR = 0.06; % min-1

p.konVA  = 5.52e-6; % pM-1 min-1
p.koffVA = 0.012;    % min-1

% compartmental transport
p.k_bt = 4.12e-4; % min-1
p.k_br = 3.19e-3; % min-1
p.k_tb = 4.12e-4; % min-1
p.k_rb = 3.19e-3+6.15e-4; % min-1; additional effect of lymphatic transport
p.AbEx = 1; % set to 0 for no antibody transport; set to 1 for regular transport

% antibody dose
dose = 10*70/150000000*1e12; % 10mg/kg * person size = 70 kg => 700 mg; MW = 150 kDa = 150,000,000 mg/mol

% simulations - 1 - run to steady state, no antibody (observe levels of
% VEGF in each tissue)

% initial conditions 
y0 = zeros(15,1);
sstime = 60*24*10;
options = odeset('MaxStep', 5e-1, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-2);
[T1,Y1] = ode15s(@VEGFAbeqns,[0:1:sstime],y0,options,p,n);
y0ss = Y1(end,:)';

p.AbEx = 0; % set to 0 for no antibody transport; set to 1 for regular transport

y0ss(n.Ab_b)=y0ss(n.Ab_b)+dose/p.Vol_b; % add antibody dose

endtime = 60*24*21; 
[T2,Y2] = ode15s(@VEGFAbeqns,[sstime:1:(sstime+endtime)],y0ss,options,p,n);
Tout = [T1;T2]-sstime; % normalize to antibody addition at time zero
Yout = [Y1;Y2];

%% VISUALIZATION

%%  Question 1
% figure('visible',on); % only show figures if flag is 'on'
f1 = figure('Name','Fig 1A');

plot( Tout/(60*24),Yout(:,n.VEGF_b),'k--','linewidth',3, 'DisplayName', 'Blood')
hold on
title( 'Reproduction of Figure 1A')
ylabel( '[VEGF] (pM)')
xlabel( 'time (days)')

legend

plot( Tout/(60*24),Yout(:,n.VEGF_t),'k:','linewidth',3,'DisplayName', 'Tumor')
% title('Concentration of free VEGF in Tumor')
ylabel( '[VEGF] (pM)')
xlabel( 'time (days)')
legend

plot( Tout/(60*24),Yout(:,n.VEGF_r),'k','linewidth',3, 'DisplayName', 'Tissue')
% title('Concentration of free VEGF in Rest of Body')
ylabel( '[VEGF] (pM)')
xlabel( 'time (days)')
legend
xlim([0 25])
hold off

f2 = figure('Name','Figure S1A');

plot(Tout/(60*24),Yout(:,n.Ab_b)*10^-6,'k--','linewidth',3, 'DisplayName', 'Blood')
hold on
title('Reproduction of Figure S1A')
ylabel('[Ab] (uM)')
xlabel('time (days)')

plot(Tout/(60*24),Yout(:,n.Ab_t)*10^-6,'k:','linewidth',3, 'DisplayName', 'Tumor')
% title('Concentration of antibody in Tumor')
ylabel('[Ab] (uM)')
xlabel('time (days)')

plot(Tout/(60*24),Yout(:,n.Ab_r)*10^-6,'k','linewidth',3, 'DisplayName', 'Tissue')
% title('Concentration of antibody in Rest of Body')
ylabel('[Ab] (uM)')
xlabel('time (days)')
xlim([0 25])
legend
hold off

f3 = figure('Name','Figure S2A');

plot(Tout/(60*24),Yout(:,n.VEGFAb_b)*10^-3,'k--','linewidth',3,'DisplayName', 'Blood')
hold on
title('Reproduction of Figure S2A')
ylabel('[Ab-V] (nM)')
xlabel('time (days)')

plot(Tout/(60*24),Yout(:,n.VEGFAb_t)*10^-3,'k:','linewidth',3, 'DisplayName', 'Tumor')
% title('Concentration of VEGF-antibody complex in Tumor')
ylabel('[Ab-V] (nM)')
xlabel('time (days)')

plot(Tout/(60*24),Yout(:,n.VEGFAb_r)*10^-3,'k','linewidth',3, 'DisplayName', 'Tissue')
% title('Concentration of VEGF-antibody complex in Rest of Body')
ylabel('[Ab-V] (nM)')
xlabel('time (days)')
xlim([0 25])
legend
hold off


p.AbEx = 1; % set to 0 for no antibody transport; set to 1 for regular transport

y0 = zeros(15,1);
sstime = 60*24*10;
options = odeset('MaxStep', 5e-1, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-2);
[T1,Y1] = ode15s(@VEGFAbeqns,[0:1:sstime],y0,options,p,n);
y0ss = Y1(end,:)';


y0ss(n.Ab_b)=y0ss(n.Ab_b)+dose/p.Vol_b; % add antibody dose
endtime = 60*24*21; 

[T2,Y2] = ode15s(@VEGFAbeqns,[sstime:1:(sstime+endtime)],y0ss,options,p,n);
Tout = [T1;T2]-sstime; % normalize to antibody addition at time zero
Yout = [Y1;Y2];


f4 = figure('Name','Figure 1B');

plot(Tout/(60*24),Yout(:,n.VEGF_b),'k--','linewidth',3, 'DisplayName', 'Blood')
hold on
title('Reproduction of Figure 1B')
ylabel('[VEGF] (pM)')
xlabel('time (days)')

plot(Tout/(60*24),Yout(:,n.VEGF_t),'k:','linewidth',3, 'DisplayName', 'Tumor')
% title('Concentration of free VEGF in Tumor')
ylabel('[VEGF] (pM)')
xlabel('time (days)')

plot(Tout/(60*24),Yout(:,n.VEGF_r),'k','linewidth',3, 'DisplayName', 'Tissue')
% title('Concentration of free VEGF in Rest of Body')
ylabel('[VEGF] (pM)')
xlabel('time (days)')
xlim([0 25])
legend
hold off

f5 = figure('Name','Figure S1B');

plot(Tout/(60*24),Yout(:,n.Ab_b)*10^-6,'k--','linewidth',3, 'DisplayName', 'Blood')
hold on
title('Reproduction of Figure S1B')
ylabel('[Ab] (uM)')
xlabel('time (days)')


plot(Tout/(60*24),Yout(:,n.Ab_t)*10^-6,'k:','linewidth',3,'DisplayName', 'Tumor')
% title('Concentration of antibody in Tumor')
ylabel('[Ab] (uM)')
xlabel('time (days)')



plot(Tout/(60*24),Yout(:,n.Ab_r)*10^-6,'k','linewidth',3, 'DisplayName', 'Tissue')
% title('Concentration of antibody in Rest of Body')
ylabel('[Ab] (uM)')
xlabel('time (days)')
xlim([0 25])
legend
hold off

f6 = figure('Name','Figure S2B');

plot(Tout/(60*24),Yout(:,n.VEGFAb_b)*10^-3,'k--','linewidth',3,'DisplayName', 'Blood')
hold on
title('Reproduction of Figure S2B')
ylabel('[Ab-V] (nM)')
xlabel('time (mins)')

plot(Tout/(60*24),Yout(:,n.VEGFAb_t)*10^-3,'k:','linewidth',3,'DisplayName', 'Tumor')
ylabel('[Ab-V] (nM)')
xlabel('time (mins)')

plot(Tout/(60*24),Yout(:,n.VEGFAb_r)*10^-3,'k','linewidth',3, 'DisplayName', 'Tissue')
ylabel('[Ab-V] (nM)')
xlabel('time (mins)')
xlim([0 25])
legend
hold off



%% Question 2



% antibody dose
dose = 10*70/150000000*1e12; % 10mg/kg * person size = 70 kg => 700 mg; MW = 150 kDa = 150,000,000 mg/mol

% simulations - 1 - run to steady state, no antibody (observe levels of
% VEGF in each tissue)

% initial conditions 
y0 = zeros(15,1);
sstime = 60*24*10;
options = odeset('MaxStep',5e-1, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-2);
[T1,Y1] = ode15s(@VEGFAbeqns,[0:1:sstime],y0,options,p,n);
y0ss = Y1(end,:)';


p.AbEx = 0; % set to 0 for no antibody transport; set to 1 for regular transport

y0ss(n.Ab_b)=y0ss(n.Ab_b); % add antibody dose

daily_dose = dose/(10*p.Vol_b);

t_curr = 0;
y_curr = y0ss;
T_all = [];
Y_all = [];

for i= 1:10
dose_time = i*60*24;
tspan = t_curr:1:dose_time;

[t_day_dose, y_day_dose] = ode15s(@VEGFAbeqns, tspan, y_curr, options, p, n);

y_curr = y_day_dose(end, :)';
y_curr(n.Ab_b) = y_curr(n.Ab_b) + daily_dose;


t_curr = t_curr + 60*24;

T_all = [T_all; t_day_dose];
Y_all = [Y_all; y_day_dose];


T_all = [T_all; t_day_dose(end) + 1e-6];
Y_all = [Y_all; y_curr'];
end


endtime = 60*24*21; 
Tout_temp = [T_all]; % normalize to antibody addition at time zero
Yout_temp = [Y_all];


[T2,Y2] = ode15s(@VEGFAbeqns,[t_curr:1:(t_curr+endtime)],y_curr,options,p,n);

Tout = [Tout_temp;T2]; % normalize to antibody addition at time zero
Yout = [Yout_temp;Y2];


f7 = figure('Name','Figure 1C');

plot( Tout/(60*24),Yout(:,n.VEGF_b),'k--','linewidth',3, 'DisplayName', 'Blood')
hold on
title( 'Reproduction of Figure 1C')
ylabel( '[VEGF] (pM)')
xlabel( 'time (days)')

plot( Tout/(60*24),Yout(:,n.VEGF_t),'k:','linewidth',3, 'DisplayName', 'Tumor')
% title('Concentration of free VEGF in Tumor')
ylabel( '[VEGF] (pM)')
xlabel( 'time (days)')


plot( Tout/(60*24),Yout(:,n.VEGF_r),'k','linewidth',3, 'DisplayName', 'Tissue')
% title('Concentration of free VEGF in Rest of Body')
ylabel( '[VEGF] (pM)')
xlabel( 'time (days)')
% xlim([0 25])
legend
hold off

f8 = figure('Name','Figure S1C');

plot(Tout/(60*24),Yout(:,n.Ab_b)*10^-6,'k--','linewidth',3, 'DisplayName', 'Blood')
hold on
title('Reproduction of Figure S1C')
ylabel('[Ab] (uM)')
xlabel('time (days)')

plot(Tout/(60*24),Yout(:,n.Ab_t)*10^-6,'k:','linewidth',3,'DisplayName', 'Tumor')
% title('Concentration of antibody in Tumor')
ylabel('[Ab] (uM)')
xlabel('time (days)')

plot(Tout/(60*24),Yout(:,n.Ab_r)*10^-6,'k','linewidth',3, 'DisplayName', 'Tissue')
% title('Concentration of antibody in Rest of Body')
ylabel('[Ab] (uM)')
xlabel('time (days)')
xlim([0 25])
legend
hold off

f9 = figure('Name','Figure S2C');

plot(Tout/(60*24),Yout(:,n.VEGFAb_b)*10^-3,'k--','linewidth',3,'DisplayName', 'Blood')
hold on
title('Reproduction of Figure S2C')
ylabel('[Ab-V] (nM)')
xlabel('time (mins)')

plot(Tout/(60*24),Yout(:,n.VEGFAb_t)*10^-3,'k:','linewidth',3,'DisplayName', 'Tumor')
ylabel('[Ab-V] (nM)')
xlabel('time (mins)')

plot(Tout/(60*24),Yout(:,n.VEGFAb_r)*10^-3,'k','linewidth',3, 'DisplayName', 'Tissue')
ylabel('[Ab-V] (nM)')
xlabel('time (mins)')
xlim([0 25])
legend
hold off

















% with extravation

dose = 10*70/150000000*1e12; % 10mg/kg * person size = 70 kg => 700 mg; MW = 150 kDa = 150,000,000 mg/mol

% simulations - 1 - run to steady state, no antibody (observe levels of
% VEGF in each tissue)

% initial conditions 
y0 = zeros(15,1);
sstime = 60*24*10;
options = odeset('MaxStep',5e-1, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-2);
[T1,Y1] = ode15s(@VEGFAbeqns,[0:1:sstime],y0,options,p,n);
y0ss = Y1(end,:)';


p.AbEx = 1; % set to 0 for no antibody transport; set to 1 for regular transport

y0ss(n.Ab_b)=y0ss(n.Ab_b); % add antibody dose

daily_dose = dose/(10*p.Vol_b);

t_curr = 0;
y_curr = y0ss;
T_all = [];
Y_all = [];

for i= 1:10
dose_time = i*60*24;
tspan = t_curr:1:dose_time;

[t_day_dose, y_day_dose] = ode15s(@VEGFAbeqns, tspan, y_curr, options, p, n);

y_curr = y_day_dose(end, :)';
y_curr(n.Ab_b) = y_curr(n.Ab_b) + daily_dose;


t_curr = t_curr + 60*24;

T_all = [T_all; t_day_dose];
Y_all = [Y_all; y_day_dose];


T_all = [T_all; t_day_dose(end) + 1e-6];
Y_all = [Y_all; y_curr'];
end


endtime = 60*24*21; 
Tout_temp = [T_all]; % normalize to antibody addition at time zero
Yout_temp = [Y_all];


[T2,Y2] = ode15s(@VEGFAbeqns,[t_curr:1:(t_curr+endtime)],y_curr,options,p,n);

Tout = [Tout_temp;T2]; % normalize to antibody addition at time zero
Yout = [Yout_temp;Y2];



f10 = figure('Name','Figure 1D');

plot(Tout/(60*24),Yout(:,n.VEGF_b),'k--','linewidth',3, 'DisplayName', 'Blood')
hold on
title('Reproduction of Figure 1D')
ylabel('[VEGF] (pM)')
xlabel('time (days)')

plot(Tout/(60*24),Yout(:,n.VEGF_t),'k:','linewidth',3, 'DisplayName', 'Tumor')
% title('Concentration of free VEGF in Tumor')
ylabel('[VEGF] (pM)')
xlabel('time (days)')

plot(Tout/(60*24),Yout(:,n.VEGF_r),'k','linewidth',3, 'DisplayName', 'Tissue')
% title('Concentration of free VEGF in Rest of Body')
ylabel('[VEGF] (pM)')
xlabel('time (days)')
xlim([0 25])
legend
hold off

f11 = figure('Name','Figure S1D');

plot(Tout/(60*24),Yout(:,n.Ab_b)*10^-6,'k--','linewidth',3, 'DisplayName', 'Blood')
hold on
title('Reproduction of Figure S1D')
ylabel('[Ab] (uM)')
xlabel('time (days)')


plot(Tout/(60*24),Yout(:,n.Ab_t)*10^-6,'k:','linewidth',3,'DisplayName', 'Tumor')
% title('Concentration of antibody in Tumor')
ylabel('[Ab] (uM)')
xlabel('time (days)')



plot(Tout/(60*24),Yout(:,n.Ab_r)*10^-6,'k','linewidth',3, 'DisplayName', 'Tissue')
% title('Concentration of antibody in Rest of Body')
ylabel('[Ab] (uM)')
xlabel('time (days)')
xlim([0 25])
legend
hold off

f12 = figure('Name','Figure S2D');

plot(Tout/(60*24),Yout(:,n.VEGFAb_b)*10^-3,'k--','linewidth',3,'DisplayName', 'Blood')
hold on
title('Reproduction of Figure S2D')
ylabel('[Ab-V] (nM)')
xlabel('time (mins)')

plot(Tout/(60*24),Yout(:,n.VEGFAb_t)*10^-3,'k:','linewidth',3,'DisplayName', 'Tumor')
ylabel('[Ab-V] (nM)')
xlabel('time (mins)')

plot(Tout/(60*24),Yout(:,n.VEGFAb_r)*10^-3,'k','linewidth',3, 'DisplayName', 'Tissue')
ylabel('[Ab-V] (nM)')
xlabel('time (mins)')
xlim([0 25])
legend
hold off




myFlag=1;   % Program finished