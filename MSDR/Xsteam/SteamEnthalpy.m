temp_table = [287.78, 315.56, 343.33]; % temp in deg-C
pres_inp = [41.37, 48.26, 55.16, 68.95, 86.18, 103.42]; % pressure in bar
pres_table = pres_inp*1e5; % bar to Pa

m = length(temp_table);
n = length(pres_table);

% matrix of zeros to store result
H_s_table = zeros(m, n); % Enthalpy of steam

for i=1:m
    for j=1:n
        H_s_table(i, j) = XSteam('h_pT', pres_inp(j), temp_table(i))/1e3; % in MJ/kg
    end
end

% Correction for error in the data in XSteam
H_s_table(1 ,5) = H_s_table(1 ,5) + 1.4;
H_s_table(1 ,6) = H_s_table(1 ,6) + 1.3;
