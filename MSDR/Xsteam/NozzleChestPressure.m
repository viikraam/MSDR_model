% density of steam in the nozzle chest in kg/m^3
rho_table = [11.2129, 16.0185, 20.8240, 25.6295, 30.4351, 35.2406];

% Enthalpy of steam in the nozzle chest
% in Imperial units -_-
H_imp = [800,810.2041,820.4082,830.6123,840.8163,851.0204,861.2245, ...
    871.4286,881.6327,891.8367,902.0408,912.2449,922.4490,932.6531, ...
    942.8571,953.0612,963.2653,973.4694,983.6735,993.8776,1004.0816,...
    1014.2857,1024.4898,1034.6939,1044.8980,1055.1021,1065.3061, ...
    1075.5102,1085.7143,1095.9184,1106.1225,1116.3265,1126.5306, ...
    1136.7347,1146.9388,1157.1429,1167.3469,1177.5510,1187.7551, ...
    1197.9592,1208.1633,1218.3673,1228.5714,1238.7755,1248.9796, ...
    1259.1837,1269.3878,1279.5918,1289.7959,1300];
% converting from imperial units [BTU/lbm] to [kJ/kg]
H_inp = H_imp*(1.05435/0.45359237);
H_table = H_inp/1e3; % from kJ/kg to MJ/kg

% length of each array
m = length(rho_table);
n = length(H_inp);

% matrix of zeros to store result
P_nc_table = zeros(n, m);

for i=1:m
    for j=1:n
        % Call XSteam function 'p_hrho' for each value of enthalpy
        % while holding density constant
        % Save the result in each row of the result matrix
        P_nc_table(j,i) = XSteam('p_hrho', H_inp(j), rho_table(i))*1e5; % bar to Pa
    end
end
    
