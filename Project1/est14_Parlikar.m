function x = est14_Parlikar(del_P, Tn, MAP, tau)
% EST14_Parlikar  CO estimator 14: Parlikar

% PPn = Pulse pressure in the nth cardiac cycle --> PPn = SAP - DAP
% SAPn = Systolic Blood Pressure
% DAPn = Diastolic Blood Pressure
% Tn = beat that begins at time tn and ends at time tn+1 --> TN = tn+1 - tn

% delta Pn = P(tn+1)-P (tn) --> beat-to-beat pressure change at the onset times


% Pn with line on top = average ABP over the nth cycle

% tau_n = RnCn



x = (del_P./Tn) + (MAP./tau);


%%%%