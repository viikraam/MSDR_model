% density of steam in the nozzle chest in kg/m^3
rho_rh_table = [1.60185, 3.20369, 4.80554, 6.40739, 8.00923, 9.61108];

% Enthalpy of steam in the nozzle chest
% in Imperial units -_-
H_rh_imp = [800,826.3158,852.6316,878.9474,905.2632,931.5790,957.8947, ...
    984.2105,1010.5263,1036.8421,1063.1579,1089.4737,1115.7895,1142.1053, ...
    1168.4211,1194.7369,1221.0526,1247.3684,1273.6842,1300];
% converting from imperial units [BTU/lbm] to [kJ/kg]
H_rh_inp = H_rh_imp*(1.05435/0.45359237);
H_rh_table = H_rh_inp/1e3; % from kJ/kg to MJ/kg

% length of each array
m = length(rho_table);
n = length(H_rh_inp);

% matrix of zeros to store result
P_rh_table = zeros(n, m);

for i=1:m
    for j=1:n
        % Call XSteam function 'p_hrho' for each value of enthalpy
        % while holding density constant
        % Save the result in each row of the result matrix
        P_rh_table(j,i) = XSteam('p_hrho', H_rh_inp(j), rho_rh_table(i))*1e5; % bar to Pascal
    end
end
    
