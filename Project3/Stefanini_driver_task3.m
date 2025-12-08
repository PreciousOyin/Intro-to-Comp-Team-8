function myFlag = ICMP_Team8_Proj3_Task3()

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
N = 100;


Vr_pt_samps = normrnd(p.Vol_r, p.Vol_r*0.25, [N,1]);
cl_a_pt_samps = normrnd(p.kcl_A, p.kcl_A*0.5, [N,1]);
cl_va_pt_samps = cl_a_pt_samps;

Vr_pt_samps(Vr_pt_samps < 0.1*p.Vol_r)       = 0.1*p.Vol_r;
cl_a_pt_samps(cl_a_pt_samps < 0.1*p.kcl_A)   = 0.1*p.kcl_A;
cl_va_pt_samps(cl_va_pt_samps < 0.1*p.kcl_A)   = 0.1*p.kcl_A;


T_all = cell(N,1);
Y_all = cell(N,1);
parfor i=1:N
fprintf('Progress: %d\n', i);
    
p_i = p;
p_i.Vol_r = Vr_pt_samps(i);
p_i.kcl_A = cl_a_pt_samps(i);
p_i.kcl_VA = cl_va_pt_samps(i);
% fprintf('%f %f\n', cl_va_pt_samps(i), cl_a_pt_samps(i));


% simulations - 1 - run to steady state, no antibody (observe levels of
% VEGF in each tissue)

% initial conditions 
y0 = zeros(15,1);
sstime = 60*24*10;
options = odeset('MaxStep',5e-1, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-2);
[T1,Y1] = ode15s(@VEGFAbeqns,[0:1:sstime],y0,options,p_i,n);
y0ss = Y1(end,:)';

p_i.AbEx = 1; % set to 0 for no antibody transport; set to 1 for regular transport

y0ss(n.Ab_b)=y0ss(n.Ab_b)+dose/p_i.Vol_b; % add antibody dose

endtime = 60*24*21; 
[T2,Y2] = ode15s(@VEGFAbeqns,[sstime:1:(sstime+endtime)],y0ss,options,p_i,n);
Tout = [T1;T2]-sstime; % normalize to antibody addition at time zero
Yout = [Y1;Y2];


T_all{i} = Tout;
Y_all{i} = Yout;

end



% VEGF and Antibody in tumor
figure;
subplot(2,2,1)
hold on;
for i = 1:N
    plot(T_all{i}/(60*24),Y_all{i}(:,n.VEGFAb_t)*10^-3,'linewidth',2)

end

title( 'Concentration of VEGF-antibody Complex in Tumor')
ylabel('[Ab-V] (nM)')
xlabel('time (days)')

xlim([0,25])
hold off


%VEGRF2 in tumor
subplot(2,2,2)
hold on;

for i = 1:N
    plot(T_all{i}/(60*24),Y_all{i}(:,n.VEGFR2_t)*10^-3,'linewidth',2)

end

title( 'Concentration of VEGF2 in Tumor')
ylabel('[[VEGFR2] (nM)')
xlabel('time (days)')

xlim([0,25])
hold off


% VEGF-VEGFR2 in Tumor
subplot(2,2,3)
hold on;

for i = 1:N
    plot(T_all{i}/(60*24),Y_all{i}(:,n.VEGFVEGFR2_t),'linewidth',2)

end

title( 'Concentration of VEGF-VEGFR2 in Tumor')
ylabel('[VEGF-VEGFR2] (pM)')
xlabel('time (days)')

xlim([0,25])
hold off


% VEGF in Tumor

subplot(2,2,4)
hold on;

for i = 1:N
    plot(T_all{i}/(60*24),Y_all{i}(:,n.VEGF_t),'linewidth',2)
end

title( 'Concentration of VEGF in Tumor')
ylabel('[VEGF] (pM)')
xlabel('time (days)')

xlim([0,25])
hold off

myFlag=1;   % Program finished