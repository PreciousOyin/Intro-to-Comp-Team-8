
function myFlag = ICMP_Team8_Proj2()

%% Question 2
clear all
close all
myFlag = 0;


ax = 0.04; % Units = min^-1
ag = 0.03; % Units = min^-1
I2C = 0.1; 
Kx = [10, 20]; % Units = U^-1


% y(1) = G(t)
% y(2) = Ug(t)
% y(3) = .Ug(t)

% y(4) = X(t)

% y(5) = X1(t)


% off set the time to account for meal duration


% insulin with the meal
Ucho = @(t) (t >= 30 && t <= 60) * (20/30);

t_bolus_start = 30;
t_bolus_end   = 60;

Ui = @(t) (t >= t_bolus_start && t <= t_bolus_end) * 2/30; 


tspan = [0 1440];
IC = [90; 0; 0; 0; 0;];

[t, y] = ode45(@(t,y) odefun(t, y, ax, ag, I2C, Kx(1), Ucho, Ui), tspan, IC);

figure
subplot(2,1,1)
hold on
plot(t, y(:,1), 'DisplayName', 'Kx = 10')
[t, y] = ode45(@(t,y) odefun(t, y, ax, ag, I2C, Kx(2), Ucho, Ui), tspan, IC);
plot(t, y(:,1), 'DisplayName', 'Kx = 20')
xlim([0 300])
ylim([84 90])
xlabel('Time [Minutes]')
ylabel('BG concentration (mg/dl)')
title('BG concentration with meal and bolus')
grid on
legend;



subplot(2,1,2)
% insulin bolus only
Ucho = @(t) (t >= 30 && t <= 60) * (0/30);
[t, y] = ode45(@(t,y) odefun(t, y, ax, ag, I2C, Kx(1), Ucho, Ui), tspan, IC);
hold on
plot(t, y(:,1), 'DisplayName', 'Kx = 10')

[t, y] = ode45(@(t,y) odefun(t, y, ax, ag, I2C, Kx(2), Ucho, Ui), tspan, IC);
plot(t, y(:,1), 'DisplayName', 'Kx = 20')
xlim([0 300])
ylim([0 95])
legend;
xlabel('Time [Minutes]')
ylabel('BG concentration (mg/dl)')
title('BG concentration with only insulin bolus')
grid on
hold off








% insulin before the meal
Ucho = @(t) (t >= 30 && t <= 60) * (20/30);

t_bolus_start = 0;
t_bolus_end   = 30;

Ui = @(t) (t >= t_bolus_start && t <= t_bolus_end) * 2/30; 


tspan = [0 1440];
IC = [90; 0; 0; 0; 0;];


[t, y] = ode45(@(t,y) odefun(t, y, ax, ag, I2C, Kx(1), Ucho, Ui), tspan, IC);

figure
subplot(2,1,1)
hold on
plot(t, y(:,1), 'DisplayName', 'Kx = 10')
[t, y] = ode45(@(t,y) odefun(t, y, ax, ag, I2C, Kx(2), Ucho, Ui), tspan, IC);
plot(t, y(:,1), 'DisplayName', 'Kx = 20')
xlim([0 300])
% ylim([84 90])
xlabel('Time [Minutes]')
ylabel('BG concentration (mg/dl)')
title('BG concentration with meal (insulin was taken before meal)')
grid on
legend;
hold off

% insulin after the meal
Ucho = @(t) (t >= 30 && t <= 60) * (20/30);

t_bolus_start = 60;
t_bolus_end   = 90;

Ui = @(t) (t >= t_bolus_start && t <= t_bolus_end) * 2/30; 


tspan = [0 1440];
IC = [90; 0; 0; 0; 0;];


[t, y] = ode45(@(t,y) odefun(t, y, ax, ag, I2C, Kx(1), Ucho, Ui), tspan, IC);

subplot(2,1,2)
hold on
plot(t, y(:,1), 'DisplayName', 'Kx = 10')
[t, y] = ode45(@(t,y) odefun(t, y, ax, ag, I2C, Kx(2), Ucho, Ui), tspan, IC);
plot(t, y(:,1), 'DisplayName', 'Kx = 20')
xlim([0 300])
% ylim([84 90])
xlabel('Time [Minutes]')
ylabel('BG concentration (mg/dl)')
title('BG concentration with meal (insulin was taken after meal)')
grid on
legend;
hold off






%% Question 3


ax = 0.2; % Units = min^-1
ag = 0.19; % Units = min^-1


% insulin only
Ucho = @(t) (t >= 0 && t <= 30) * (0);

t_bolus_start = 30;
t_bolus_end   = 60;

Ui = @(t) (t >= t_bolus_start && t <= t_bolus_end) * (2/30); 


tspan = [0 1440];
IC = [90; 0; 0; 0; 0;];


% [t, y] = ode45(@(t,y) odefun(t, y, ax, ag, I2C, Kx(1), Ucho, Ui), tspan, IC);

figure
subplot(2,1,1)
hold on
% plot(t, y(:,1), 'DisplayName', 'Kx = 10')
[t, y] = ode45(@(t,y) odefun(t, y, ax, ag, I2C, Kx(2), Ucho, Ui), tspan, IC);
plot(t, y(:,1), 'DisplayName', 'Kx = 20', 'LineWidth', 2)
xlim([0 300])
% ylim([84 90])
xlabel('Time [Minutes]')
ylabel('BG concentration (mg/dl)')
title('BG concentration (insulin bolus only)')
grid on
% legend;
hold off









% meal only
Ucho = @(t) (t >= 30 && t <= 60) * (20/30);

t_bolus_start = 60;
t_bolus_end   = 90;

Ui = @(t) (t >= t_bolus_start && t <= t_bolus_end) * 0; 


tspan = [0 1440];
IC = [90; 0; 0; 0; 0;];


% [t, y] = ode45(@(t,y) odefun(t, y, ax, ag, I2C, Kx(1), Ucho, Ui), tspan, IC);
% figure
subplot(2,1,2)
hold on
% plot(t, y(:,1), 'DisplayName', 'Kx = 10')
[t, y] = ode45(@(t,y) odefun(t, y, ax, ag, I2C, Kx(2), Ucho, Ui), tspan, IC);
plot(t, y(:,1), 'DisplayName', 'Kx = 20', 'LineWidth', 2)
xlim([0 300])
% ylim([84 90])
xlabel('Time [Minutes]')
ylabel('BG concentration (mg/dl)')
title('BG concentration (meal only)')
grid on
% legend;
hold off


myFlag=1;   % Program finished



function dydt = odefun(t, y, ax, ag, I2C, Kx, Ucho, Ui)
    dydt = zeros (5,1);
    dydt(1) = -y(4) + y(2);
    dydt(2) = y(3);
    dydt(3) = (-2*ag*y(3)) - ((ag^2)*y(2)) + I2C*Kx*(ag^2)*Ucho(t);
    dydt(4) = -ax*y(4) + ax*y(5);
    dydt(5) = -ax*y(5)+ Kx*ax*Ui(t);




return