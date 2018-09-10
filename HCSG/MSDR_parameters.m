% MSDR parameters file
% Authors: Vikram Singh, Alex Wheeler
% Supervisors: Dr. Ondrej Chvala, Dr. Belle Upadhyaya
% Project: ORNL/UTK MSR dynamic modeling

%% CORE PARAMETERS
%% Core heat transfer parameters

core_height = 6.8; % (m), Source: ORNL-3832, table 9
core_radius = 3.4; % (m), Source: ORNL-3832, table 9
core_fuel_frac = 0.1; % Source: ORNL-3832, table 9

core_volume = pi*(core_radius)^2*core_height; %(m3) 
core_fuelvol = 0.1*core_volume; 
% ^^THERE IS A DIFFERENT VALUE QUOTED ABOVE, This one is calculated

% fuel salt density
den_f   = 3355.95           ; % fuel salt density as calc. using empirical eqn. 236.3-(2.33e-2*T(deg-F)), (kg/m3), Source: Source: ORNL-TM-3832, table 3
den_g   = 1860              ; % density of graphite (kg/m3) same as MSRE, Source: ORNL-TM-0728, table 5.3

mass_fuel = core_fuelvol*den_f; % (kg)
mass_fn = mass_fuel/4;
mass_graphite = (core_volume - core_fuelvol)*den_g; % (kg)
mass_gn = mass_graphite/4;


% fuel salt specific heat capacity, Source: ORNL-TM-4676, table 3.1
Cp_f    = 1.355616e-3         ; % (MJ/kg/°C)

% graphite specific heat capacity, Source: Same as MSRE
Cp_g    = 1.773e-3          ; % (MJ/kg/°C)

% nodal masses in core
mfa     = [mass_fn, mass_fn]; % fuel nodes a
mfb     = [mass_fn, mass_fn]; % fuel nodes b
mg      = [mass_gn, mass_gn]; % graphite nodes

% heat deposition fractions, adopted from MSRE
ka       = [0.235, 0.235]; % power generation fraction for fuel nodes a
kb       = [0.235, 0.235]; % power generation fraction for fuel nodes b
kga      = [0.0075, 0.0075]; % heat transferred from for graphite node to fuel node a
kgb      = [0.0075, 0.0075]; % heat transferred from for graphite node to fuel node b

% (mass)*(specific heat capacity)
mcp_fa   = mfa*Cp_f         ;
mcp_fb   = mfb*Cp_f         ;
mcp_g    = mg*Cp_g          ; 

% (heat transfer coefficient)*(area) for fuel/graphite interfaces 
h_fg     = 1.8e-4           ; % (MW/m2/°C)
A_fg     = 8979             ; % (m2)
hA_fg    = [h_fg*A_fg/2, h_fg*A_fg/2]; 

% Importance factors for core nodes
Ifa      = [0.25, 0.25]; % fuel nodes a
Ifb      = [0.25, 0.25]; % fuel nodes b
Ig       = [0.5, 0.5]; % graphite nodes

% initial temperatures (°C)
% Fuel volumes
core_fv = core_fuelvol      ; % volume of fuel in the core (m3), Source: ORNL-TM-3832, table 9
pipe_fv =  9.0614           ; % volume of fuel in external piping and PHX (m3), Source: ORNL-TM-3832, table 9
tot_fv  = 30.5822           ; % total volume of fuel salt in the plant (m3), Source: ORNL-TM-3832, table 9

% Fuel temperature in core
Tf_out  = 676.667           ; % core fuel salt outlet temperature (°C), Source: ORNL-TM-3832, fig 1 & table 2
Tf_in   = 565.556           ; % core fuel salt inlet temperature (°C), Source: ORNL-TM-3832, fig 1 & table 2
Tf_avg  = (Tf_out + Tf_in)/2; % core fuel salt avg temperature, (out+in)/2

T0_fa    = [593.3337, 648.8891]; % fuel nodes a
T0_fb    = [621.1114, Tf_out]; % fuel nodes b
T0_g     = [607.255, 662.811]; % graphite nodes

%% PKE tau calculations
P       = 750; % (MW)

% mass of fuel salt
core_fm = core_fv*den_f     ; % mass of fuel in core (kg)
pipe_fm = pipe_fv*den_f     ; % mass of fuel in external pipes and PHX (kg)

% fuel salt flow rate
W_f     = 4989.48          ; % fuel salt flow rate (kg/s)

W_f1    = W_f/3;
W_f2    = W_f/3;
W_f3    = W_f/3;

m_fm    = W_f*2; 

% Source: ORNL-TM-3832, table 2, and Source: ORNL-4676 table 3.2, lists parameters for a 150 MW(th) HX, there
% are six of these in the MSDR for a total of 750 MW(th)

% PKE taus
tau_c   = core_fm/W_f       ; % fuel salt core residence time (s)
tau_l   = pipe_fm/W_f       ; % fuel salt loop residence time (s)

tau_p = 5;
therm_conv = 0.05;


% reactdata = [0 2e-4 0 -2e-4 0 2e-4 0 -2e-4 0 2e-4 0];
% reacttime = [0 2500 2550 2600 2650 2700 2750 2800 2850 2900 2950];
reactdata = [0 (-20e-5)];
reacttime = [0 0e3];
react = timeseries(reactdata,reacttime);

%% PKE parameters

% U-235 parameters, Source: ORNL-TM-2571, table 1
Lam     = 2.400E-04         ; % mean generation time (s) 
lam     = [0.0124, 0.0305, 0.111, 0.301, 1.14, 3.01]; % group decay constants (s^-1)
beta    = [0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.000280]; % group delayed-neutron yield (abs frac)
beta_t  = sum(beta)         ; % total delayed-neutron yield (abs frac)

rho_0   = beta_t - bigterm(beta, lam, tau_l, tau_c); % reactivity needed for going from solid to fluid fuel (drho)

n0      = 1                 ; % initial fractional neutron density

% initial fractional delayed-neutron densities for groups 1..6
C0(1)   = ((beta(1))/Lam)*(1.0/(lam(1) - (exp(-lam(1)*tau_l) - 1.0)/tau_c));
C0(2)   = ((beta(2))/Lam)*(1.0/(lam(2) - (exp(-lam(2)*tau_l) - 1.0)/tau_c));
C0(3)   = ((beta(3))/Lam)*(1.0/(lam(3) - (exp(-lam(3)*tau_l) - 1.0)/tau_c));
C0(4)   = ((beta(4))/Lam)*(1.0/(lam(4) - (exp(-lam(4)*tau_l) - 1.0)/tau_c));
C0(5)   = ((beta(5))/Lam)*(1.0/(lam(5) - (exp(-lam(5)*tau_l) - 1.0)/tau_c));
C0(6)   = ((beta(6))/Lam)*(1.0/(lam(6) - (exp(-lam(6)*tau_l) - 1.0)/tau_c));

%% Reactivity feedback coefficients
% Source: ORNL-4541

a_f     = -3.22E-5          ; % total fuel salt thermal reactivity coefficient (drho/°C)
a_g     = +2.35E-5          ; % total graphite thermal reactivity coefficient (drho/°C)


%% HX PARAMETERS
%% PHX heat transfer parameters

Wp_phx1   =  831.58            ; % primary fluid mass flow rate in PHX, 1/6 of primary flow rate, table 2 ORNL-TM-3832 (kg/s)
Wp_phx2   =  831.58            ;
Wp_phx3   =  831.58            ;
Wp_phx4   =  831.58            ;
Wp_phx5   =  831.58            ;
Wp_phx6   =  831.58            ;

mp_phxm    = (Wp_phx1 + Wp_phx2)*1; % mass in the primary fluid mixing plenum (residence time of 1 sec)

Ws_phx1   =  466.19            ; % secondary fluid mass flow rate in PHX, table 4 ORNL-TM-3832 (kg/s)
Ws_phx2   =  466.19            ;
Ws_phx3   =  466.19            ;
Ws_phx4   =  466.19            ;
Ws_phx5   =  466.19            ;
Ws_phx6   =  466.19            ;

ms_phxm    = (Ws_phx1 + Ws_phx2)*1; % mass in the secondary fluid mixing plenum (residence time of 1 sec)

Cp_p_phx =  1339.78e-6            ; % PHX primary fluid specific heat capacity, table 2 ORNL-TM-3832 (MJ/kg/°C) 
Cp_s_phx =  2386.48e-6            ; % PHX secondary fluid specific heat capacity, table 2 ORNL-TM-3832(MJ/kg/°C)
Cp_t_phx =  577.78e-6             ; % PHX Hastelloy tube specific heat capacity, ORNL-TM-0728 pp. 20(MJ/kg/°C)

% mass of PHX nodes (kg)
% All masses were calculated using the volume in ORNL-TM-3832 table 2 and the density calculated from ORNL-TM-3832 table 3
% The mass of the tube was found using the volume in ORNL-TM-3832 table 2 and a density of 8774.5 kg/m^3 from ORNL-TM-0728 pp. 20

mp_phx1   = [494.7, 494.7, 494.7, 494.7]; % PHX1 primary nodes (kg)
mt_phx1   = [1323.6, 1323.6]; % PHX1 tube nodes (kg)
ms_phx1   = [1152.5, 1152.5, 1152.5, 1152.5]; % PHX1 secondary nodes (kg)

mp_phx2   = [494.7, 494.7, 494.7, 494.7]; % PHX2 primary nodes (kg)
mt_phx2   = [1323.6, 1323.6]; % PHX2 tube nodes (kg)
ms_phx2   = [1152.5, 1152.5, 1152.5, 1152.5]; % PHX2 secondary nodes (kg)

mp_phx3   = [494.7, 494.7, 494.7, 494.7]; % PHX3 primary nodes (kg)
mt_phx3   = [1323.6, 1323.6]; % PHX3 tube nodes (kg)
ms_phx3   = [1152.5, 1152.5, 1152.5, 1152.5]; % PHX3 secondary nodes (kg)

mp_phx4   = [494.7, 494.7, 494.7, 494.7]; % PHX4 primary nodes (kg)
mt_phx4   = [1323.6, 1323.6]; % PHX4 tube nodes (kg)
ms_phx4   = [1152.5, 1152.5, 1152.5, 1152.5]; % PHX4 secondary nodes (kg)

mp_phx5   = [494.7, 494.7, 494.7, 494.7]; % PHX5 primary nodes (kg)
mt_phx5   = [1323.6, 1323.6]; % PHX5 tube nodes (kg)
ms_phx5   = [1152.5, 1152.5, 1152.5, 1152.5]; % PHX5 secondary nodes (kg)

mp_phx6   = [494.7, 494.7, 494.7, 494.7]; % PHX6 primary nodes (kg)
mt_phx6   = [1323.6, 1323.6]; % PHX6 tube nodes (kg)
ms_phx6   = [1152.5, 1152.5, 1152.5, 1152.5]; % PHX6 secondary nodes (kg)

% (mass)*(specific heat capacity)
mcp_p_phx1 = Cp_p_phx*mp_phx1; % PHX1 primary nodes
mcp_t_phx1 = Cp_t_phx*mt_phx1; % PHX1 tube nodes
mcp_s_phx1 = Cp_s_phx*ms_phx1; % PHX1 secondary nodes

mcp_p_phx2 = Cp_p_phx*mp_phx2; % PHX2 primary nodes
mcp_t_phx2 = Cp_t_phx*mt_phx2; % PHX2 tube nodes
mcp_s_phx2 = Cp_s_phx*ms_phx2; % PHX2 secondary nodes

mcp_p_phx3 = Cp_p_phx*mp_phx3; % PHX3 primary nodes
mcp_t_phx3 = Cp_t_phx*mt_phx3; % PHX3 tube nodes
mcp_s_phx3 = Cp_s_phx*ms_phx3; % PHX3 secondary nodes

mcp_p_phx4 = Cp_p_phx*mp_phx4; % PHX4 primary nodes
mcp_t_phx4 = Cp_t_phx*mt_phx4; % PHX4 tube nodes
mcp_s_phx4 = Cp_s_phx*ms_phx4; % PHX4 secondary nodes

mcp_p_phx5 = Cp_p_phx*mp_phx5; % PHX5 primary nodes
mcp_t_phx5 = Cp_t_phx*mt_phx5; % PHX5 tube nodes
mcp_s_phx5 = Cp_s_phx*ms_phx5; % PHX5 secondary nodes

mcp_p_phx6 = Cp_p_phx*mp_phx6; % PHX6 primary nodes
mcp_t_phx6 = Cp_t_phx*mt_phx6; % PHX6 tube nodes
mcp_s_phx6 = Cp_s_phx*ms_phx6; % PHX6 secondary nodes

% heat transfer coefficient in PHX
% overall heat transfer coefficient was given, for now assumed to be the same heat transfer coefficient
h_p_phx = 7944.26e-6; % heat transfer coefficient of the primary nodes ORNL-TM-3832 table 2 (MW/m^2/°C)
h_s_phx = 7944.26e-6; % heat transfer coefficient of the secondary nodes ORNL-TM-3832 table 2 (MW/m^2/°C)

% area of heat transfer in PHX
A_p_phx = [76.81, 76.81, 76.81, 76.81]; % area of heat transfer for primary nodes ORNL-TM-3832 table 2 (m^2)
A_s_phx = [93.57, 93.57, 93.57, 93.57]; % area of heat transfer for secondary nodes ORNL-TM-3832 table 2 (m^2)

% (heat transfer coefficient)*(area)
hA_p_phx1 = h_p_phx*A_p_phx; % for primary/tube interface
hA_p_phx2 = h_p_phx*A_p_phx;
hA_p_phx3 = h_p_phx*A_p_phx;
hA_p_phx4 = h_p_phx*A_p_phx;
hA_p_phx5 = h_p_phx*A_p_phx;
hA_p_phx6 = h_p_phx*A_p_phx;

hA_s_phx1 = h_s_phx*A_s_phx; % for tube/secondary interface
hA_s_phx2 = h_s_phx*A_s_phx;
hA_s_phx3 = h_s_phx*A_s_phx;
hA_s_phx4 = h_s_phx*A_s_phx;
hA_s_phx5 = h_s_phx*A_s_phx;
hA_s_phx6 = h_s_phx*A_s_phx;

% initial temperatures (°C)
% calculated from temperatures in table 2 in ORNL-TM-3832

T0_p_phx1 = [648.89, 621.12, 593.34, 565.56]; % PHX1 primary nodes
T0_t_phx1 = [603.12, 547.57]; % PHX1 tube nodes
T0_s_phx1 = [510, 537.78, 565.55, 593.33]; % PHX1 secondary nodes

T0_p_phx2 = [648.89, 621.12, 593.34, 565.56]; % PHX2 primary nodes
T0_t_phx2 = [603.12, 547.57]; % PHX2 tube nodes
T0_s_phx2 = [510, 537.78, 565.55, 593.33]; % PHX2 secondary nodes

T0_p_phx3 = [648.89, 621.12, 593.34, 565.56]; % PHX3 primary nodes
T0_t_phx3 = [603.12, 547.57]; % PHX3 tube nodes
T0_s_phx3 = [510, 537.78, 565.55, 593.33]; % PHX3 secondary nodes

T0_p_phx4 = [648.89, 621.12, 593.34, 565.56]; % PHX4 primary nodes
T0_t_phx4 = [603.12, 547.57]; % PHX4 tube nodes
T0_s_phx4 = [510, 537.78, 565.55, 593.33]; % PHX4 secondary nodes

T0_p_phx5 = [648.89, 621.12, 593.34, 565.56]; % PHX5 primary nodes
T0_t_phx5 = [603.12, 547.57]; % PHX5 tube nodes
T0_s_phx5 = [510, 537.78, 565.55, 593.33]; % PHX5 secondary nodes

T0_p_phx6 = [648.89, 621.12, 593.34, 565.56]; % PHX6 primary nodes
T0_t_phx6 = [603.12, 547.57]; % PHX6 tube nodes
T0_s_phx6 = [510, 537.78, 565.55, 593.33]; % PHX6 secondary nodes

tau_c_hx  = .3575;
tau_hx_c  = .3575;

tau_phx_shx = 1;
%% SHX heat transfer parameters

Wp_shx1   = Ws_phx1           ; % primary fluid mass flow rate in SHX, same as PHX secondary flow rate
Wp_shx2   = Ws_phx2           ;
Wp_shx3   = Ws_phx3           ;
Wp_shx4   = Ws_phx4           ;
Wp_shx5   = Ws_phx5           ;
Wp_shx6   = Ws_phx6           ;

mp_shxm    = (Wp_shx1 + Wp_shx2)*1;  % mass in the primary fluid mixing plenum (residence time of 1 sec)

Ws_shx1   = 7.181839439165465e+02            ; % secondary fluid mass flow rate in SHX, ORNL-TM-3832 table 4 (kg/s)
Ws_shx2   = 7.181839439165465e+02            ;
Ws_shx3   = 7.181839439165465e+02            ;
Ws_shx4   = 7.181839439165465e+02            ;
Ws_shx5   = 7.181839439165465e+02            ;
Ws_shx6   = 7.181839439165465e+02            ;

W_3       = Ws_shx1 + Ws_shx2 + Ws_shx3 + Ws_shx4 + Ws_shx5 + Ws_shx6; % Hitec salt total flow rate

% Hitec salt flow rate from each loop
W_31      = W_3/3; 
W_32      = W_3/3;
W_33      = W_3/3;

m_3m      = W_3*2; % Hitec salt mixing plenum (residence time of 2 second) 

ms_shxm   = (Ws_shx1 + Ws_shx2)*1;  % mass in the secondary fluid mixing plenum (residence time of 1 sec)

Cp_p_shx = Cp_s_phx           ; % PHX primary fluid specific heat capacity (J/kg/°C) 
Cp_s_shx = 1549.12e-6            ; % PHX secondary fluid specific heat capacity ORNL-TM-3832 table 3(MJ/kg/°C)

% mass of PHX nodes (kg)
% All masses were calculated using the volume in ORNL-TM-3832 table 4 and the density calculated from ORNL-TM-3832 table 3
% The mass of the tube was found using the volume in ORNL-TM-3832 table 4 and a density of 8774.5 kg/m^3 from ORNL-TM-0728 pp. 20

mp_shx1   = [429.3, 429.3, 429.3, 429.3]; % SHX1 primary nodes
mt_shx1   = [1940, 1940]; % SHX1 tube nodes
ms_shx1   = [1861.5, 1861.5, 1861.5, 1861.5]; % SHX1 secondary nodes

mp_shx2   = [429.3, 429.3, 429.3, 429.3]; % SHX2 primary nodes
mt_shx2   = [1940, 1940]; % SHX2 tube nodes
ms_shx2   = [1861.5, 1861.5, 1861.5, 1861.5]; % SHX2 secondary nodes

mp_shx3   = [429.3, 429.3, 429.3, 429.3]; % SHX3 primary nodes
mt_shx3   = [1940, 1940]; % SHX3 tube nodes
ms_shx3   = [1861.5, 1861.5, 1861.5, 1861.5]; % SHX3 secondary nodes

mp_shx4   = [429.3, 429.3, 429.3, 429.3]; % SHX4 primary nodes
mt_shx4   = [1940, 1940]; % SHX4 tube nodes
ms_shx4   = [1861.5, 1861.5, 1861.5, 1861.5]; % SHX4 secondary nodes

mp_shx5   = [429.3, 429.3, 429.3, 429.3]; % SHX5 primary nodes
mt_shx5   = [1940, 1940]; % SHX5 tube nodes
ms_shx5   = [1861.5, 1861.5, 1861.5, 1861.5]; % SHX5 secondary nodes

mp_shx6   = [429.3, 429.3, 429.3, 429.3]; % SHX6 primary nodes
mt_shx6   = [1940, 1940]; % SHX6 tube nodes
ms_shx6   = [1861.5, 1861.5, 1861.5, 1861.5]; % SHX6 secondary nodes

% (mass)*(specific heat capacity)
mcp_p_shx1 = Cp_p_shx*mp_shx1; % SHX1 primary nodes
mcp_t_shx1 = Cp_t_phx*mt_shx1; % SHX1 tube nodes
mcp_s_shx1 = Cp_s_shx*ms_shx1; % SHX1 secondary nodes

mcp_p_shx2 = Cp_p_shx*mp_shx2; % SHX2 primary nodes
mcp_t_shx2 = Cp_t_phx*mt_shx2; % SHX2 tube nodes
mcp_s_shx2 = Cp_s_shx*ms_shx1; % SHX2 secondary nodes

mcp_p_shx3 = Cp_p_shx*mp_shx3; % SHX3 primary nodes
mcp_t_shx3 = Cp_t_phx*mt_shx3; % SHX3 tube nodes
mcp_s_shx3 = Cp_s_shx*ms_shx1; % SHX3 secondary nodes

mcp_p_shx4 = Cp_p_shx*mp_shx4; % SHX4 primary nodes
mcp_t_shx4 = Cp_t_phx*mt_shx4; % SHX4 tube nodes
mcp_s_shx4 = Cp_s_shx*ms_shx1; % SHX4 secondary nodes

mcp_p_shx5 = Cp_p_shx*mp_shx5; % SHX5 primary nodes
mcp_t_shx5 = Cp_t_phx*mt_shx5; % SHX5 tube nodes
mcp_s_shx5 = Cp_s_shx*ms_shx1; % SHX5 secondary nodes

mcp_p_shx6 = Cp_p_shx*mp_shx6; % SHX6 primary nodes
mcp_t_shx6 = Cp_t_phx*mt_shx6; % SHX6 tube nodes
mcp_s_shx6 = Cp_s_shx*ms_shx1; % SHX6 secondary nodes

% heat transfer coefficient in PHX
% overall heat transfer coefficient was given, for now assumed to be the same heat transfer coefficient
h_p_shx = 5674.46e-6; % heat transfer coefficient of the primary nodes ORNL-TM-3832 table 4 (MW/m^2/°C)
h_s_shx = 5674.46e-6; % heat transfer coefficient of the secondary nodes ORNL-TM-3832 table 4 (MW/m^2/°C)

% area of heat transfer in PHX
A_p_shx = [112.28, 112.28, 112.28, 112.28]; % area of heat transfer for primary nodes ORNL-TM-3832 table 4 (m^2)
A_s_shx = [137.142, 137.142, 137.142, 137.142]; % area of heat transfer for secondary nodes ORNL-TM-3832 table 4 (m^2)

% (heat transfer coefficient)*(area)
hA_p_shx1 = h_p_shx*A_p_shx; % for primary/tube interface
hA_p_shx2 = h_p_shx*A_p_shx;
hA_p_shx3 = h_p_shx*A_p_shx;
hA_p_shx4 = h_p_shx*A_p_shx;
hA_p_shx5 = h_p_shx*A_p_shx;
hA_p_shx6 = h_p_shx*A_p_shx;

hA_s_shx1 = h_s_shx*A_s_shx; % for tube/secondary interface
hA_s_shx2 = h_s_shx*A_s_shx;
hA_s_shx3 = h_s_shx*A_s_shx;
hA_s_shx4 = h_s_shx*A_s_shx;
hA_s_shx5 = h_s_shx*A_s_shx;
hA_s_shx6 = h_s_shx*A_s_shx;

% initial temperatures (°C)
T0_p_shx1 = [565.55, 537.78, 510, 482.22]; % SHX1 primary nodes
T0_t_shx1 = [538, 519.7]; % SHX1 tube nodes
T0_s_shx1 = [454.45, 482.23, 510, 537.78]; % SHX1 secondary nodes

T0_p_shx2 = [565.55, 537.78, 510, 482.22]; % SHX2 primary nodes
T0_t_shx2 = [538, 519.7]; % SHX2 tube nodes
T0_s_shx2 = [454.45, 482.23, 510, 537.78]; % SHX2 secondary nodes

T0_p_shx3 = [565.55, 537.78, 510, 482.22]; % SHX3 primary nodes
T0_t_shx3 = [538, 519.7]; % SHX3 tube nodes
T0_s_shx3 = [454.45, 482.23, 510, 537.78]; % SHX3 secondary nodes

T0_p_shx4 = [565.55, 537.78, 510, 482.22]; % SHX4 primary nodes
T0_t_shx4 = [538, 519.7]; % SHX4 tube nodes
T0_s_shx4 = [454.45, 482.23, 510, 537.78]; % SHX4 secondary nodes

T0_p_shx5 = [565.55, 537.78, 510, 482.22]; % SHX5 primary nodes
T0_t_shx5 = [538, 519.7]; % SHX5 tube nodes
T0_s_shx5 = [454.45, 482.23, 510, 537.78]; % SHX5 secondary nodes

T0_p_shx6 = [565.55, 537.78, 510, 482.22]; % SHX6 primary nodes
T0_t_shx6 = [538, 519.7]; % SHX6 tube nodes
T0_s_shx6 = [454.45, 482.23, 510, 537.78]; % SHX6 secondary nodes

tau_shx_phx = 1;
tau_shx_thx = 1;

T3_out       = 537.778;

%% Boiler and Reheater Parameters
mn_ms = ms_shxm; % mass salt in mixing node (kg)
p_b = 0.87*P/3; % power removed by boiler (MW)
p_r = 0.13*P/3; % power removed by reheater (MW)

W_sb = (Ws_shx1 + Ws_shx2); %p_b/(111.11*Cp_s_shx); % Coolant salt flow through boiler (kg/s)
m_sb = 1717; %4808;
W_sr = p_r/(111.11*Cp_s_shx); % Coolant salt flow through reheater (kg/s)
m_sr = 3062;

mb_cpb = 2.689619410657651e+00; %3.297E+01;
mr_cpr = 1.050E+01;

tau_fehx_b = 13.5; % (sec) fertile hx to boiler
tau_fehx_r = 17.3; % (sec) fertile hx to reheater
tau_b_fhx = 4.2; % (sec) boiler to fuel hx
tau_r_fhx = 11.1; % (sec) reheater to fuel hx

% %% TEST CONSTANT HEAT REMOVAL
% 
% rho_rs     = 810.5; % water density @255-C and 150 bar (kg/m^3)
% W_rs       = 4.86351819; % water flow rate (kg/s) 
% 
% mn_rs       = 7446; % coolant air mass in rad (kg) 
% mn_rp      = 1717; % hitec in rad (kg)
% scp_rs     = 4.92E-3; % (MJ/kg-C) specific heat capacity of air at (air_out+air_in)/2 REFPROP
% 
% A_rad    = 480; % (m^2) surface area of radiator ORNL-TM-0728 p.14
% h_roverall = P/A_rad/((800+537.78)/2-(482.222+254.5)/2); % (MW/m^2-C) polimi thesis
% 
% % Primary Side 
% mcp_rpn   = mn_rp*Cp_s_shx; % (mass of material x heat capacity of material) of fuel salt per lump in MW-s/°C
% hA_rpn    = h_roverall*A_rad; %3.030; % (primary to secondary heat transfer coeff x heat transfer area) in MW/°C
% 
% % Secondary Side - DONE
% mcp_rsn   = mn_rs*scp_rs; % (mass of material x heat capacity of material) of coolant salt per lump in MW-s/°C
% hA_rsn = h_roverall*A_rad; % (tube to secondary heat transfer coeff x heat transfer area) in MW/°C
% 

%% Xenon reactivity effects

gamma_I  = 5.1135e-2; % weighted yield for I-135
gamma_Xe = 1.1628e-2; % weighted yield for Xe-135

lam_I    = 2.875e-5;  % decay constant for I-135 (s^-1)
lam_Xe   = 2.0916e-5; % decay constant for Xe-135 (s^-1)

lam_bubl = 2.0e-2;    % effective bubbling out constant (s^-1)

sig_Xe   = 2.66449e-18; % (cm^2) microscopic cross-section for Xe (n,gamma) reaction 

molc_wt  = .715*(7.016+18.998)+.16*(9.012+2*18.998)+.12*(4*18.998+232.038)+.005*(4*18.998+235.044); % (g/mol)
molc_den = 0.001*den_f/molc_wt;          % (mol/cm^3)
U_den    = .005*molc_den*6.022E23;       % (#U/cm^3)
U_sig    = 5.851e-22;                    % (cm^2)

Sig_f_msdr = U_den*U_sig;                % (cm^-1)

phi_0 = P/(3.04414e-17*1e6*Sig_f_msdr);  % neutrons cm^-2 s^-1

I_0 = gamma_I*Sig_f_msdr*phi_0/lam_I; 
Xe_0 = (gamma_Xe*Sig_f_msdr*phi_0 + lam_I*I_0)/(lam_Xe + sig_Xe*phi_0 + lam_bubl);

Xe_og = lam_bubl*Xe_0/(lam_Xe); % initial Xe conc. in off-gas system

Sig_a = 1.02345; % (cm^-1) macroscopic absorption cross-section for core

rhoXe_0 = (-sig_Xe/Sig_a)*(gamma_I+gamma_Xe)*Sig_f_msdr*phi_0/(lam_Xe + sig_Xe*phi_0 + lam_bubl);


%% heat transfer coefficient lookup table

xdat=[0.00483288502733
0.024472562101481
0.044096021440641
0.058773071607534
0.083294286913988
0.102852875313184
0.127341655149656
0.146900243548852
0.161577293715744
0.18113588211494
0.19578049681185
0.215339085211046
0.229999917642947
0.254472479744428
0.274031068143625
0.293589656542821
0.313132027207026
0.337588371573515
0.371824010139604
0.406075866440683
0.44029528727178
0.47456336130785
0.513745409046206
0.543083291645
0.58717931308564
0.626328925354014
0.660629434860066
0.694897508896137
0.734079556634493
0.773229168902867
0.822239164045793
0.856507238081864
0.900587041787514
0.949564601460459
0.983865110966511
1.02306337643986
1.05736388594591
1.08675042174968
1.12598112269301
1.16029784993405
1.21422181877926
1.26326424939217
1.3024949503355
1.3466558427161
1.39571449106401
1.43497762747732
1.50863857980663];

ydat=[0.013508209475996
0.024765050705993
0.038273260181989
0.054032837903984
0.072043783871979
0.094557466331972
0.117071148791966
0.139584831251959
0.155344408973954
0.177858091433948
0.198120405647942
0.220634088107936
0.23864503407593
0.263410084781923
0.285923767241916
0.308437449701909
0.333202500407902
0.360218919359894
0.398492179541882
0.434514071477873
0.47503869990586
0.50880922359585
0.544831115531841
0.578601639221831
0.616874899403819
0.657399527831807
0.686667315029798
0.720437838719788
0.756459730655778
0.796984359083766
0.837508987511754
0.871279511201745
0.911804139629733
0.95683150454972
0.98609929174771
1.0198698154377
1.04913760263569
1.07615402158768
1.10542180878568
1.13243822773767
1.17521422441165
1.21123611634765
1.24050390354564
1.26977169074363
1.30354221443362
1.32830726513961
1.3688318935676];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Steam Generator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATERIAL & FLUID PROPERTIES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%from mm to ft;
%c1=0.1*0.3937/12;
c1 = 1e-3; % mm to m

% Number of tubes=820;
N = 820; 
% Tube inner diameter=13.24mm
Ri_thermal = 13.24/2*c1; % [m]
% Tube outside diameter=17.46mm
Ro_thermal = 17.46/2*c1; % [m]
% Inlet temperature on the primary side
T_pin = 328.4; %[C]

% Tube inside flow area;
Ri_hydraulic = Ri_thermal;
A_s = pi*Ri_hydraulic^2;

% Internal shell external diameter Di=0.61 m;
Di = 0.61;
% External shell internal diameter Dt=1.62m
Do = 1.62;
% Radial pitch=25mm;
t = 25*c1;
% Total shell side projected area
Ap_total = 1/4*pi*(Do^2 - Di^2)*(1 - Ro_thermal*2/t);
% Shell side hydraulic radius=2*flow area/wetting perimeter
% Wetting area
Peri_wet = N*2*pi*Ro_thermal;
Ro_hydraulic = 2*Ap_total/Peri_wet;
A_p = pi*Ro_hydraulic^2;

% Tube side hydraulic dimater
Ri_hydraulic = Ri_thermal;
% Flow area on the secondary side;
A_s = pi*(Ri_hydraulic^2);  
% Cross section for the tube;
A_w = pi*(Ro_thermal^2 - Ri_thermal^2);

% Specific heat capacity;
Cp_p        = 5591.471400012e-6; % [MJ/(kg°C)]             
Cp_w        = 456.361200001e-6;  % [MJ/(kg°C)] inconel                        
Cp_fw       = 4697.58960001e-6;  % [[MJ/(kg/°C)] feedwater                      
Cp_sc       = 4697.58960001e-6;  % [[MJ/(kg/°C)] subcooled fluid                
% Cp_s        = 3188.208e-6; % [[MJ/(kg/°C)] steam
rho_p       = 716.826235984;     % Water [kg/m^3] calculated using 130.6-(2.54e-2*T) (where T is in °F) Table3 ORNL-TM-3832  44.75 primary fluid                  
rho_w       = 8425.711734699;    % [kg/m^3] inconel                        
rho_fw      = 828.314741067;     % [kg/m^3] feedwater                      
rho_f       = 751.426116872;     % [kg/m^3] subcooled boiling density          
rho_sc      = (rho_fw + rho_f)/2;  % [kg/m^3] subcooled fluid density        
     
P_s   = 5.8; % [MPa]
T_s   = (280 + 318)/2 + 273; 
[dum,dum,Cp_s,vss] = hsh(P_s, T_s); % vss == specific volume
Cp_s  = Cp_s*1e-6; % [MJ/kg-C]
rho_s = 1/vss; % [kg/m^3]

%%%%%%%%%%%%%%%%%%%%
%%% SG TUBE SIDE %%%
%%%%%%%%%%%%%%%%%%%%

% Heat transfer coefficient
h_pw  = 9101.343578935e-6; % [MW/m^2-C]
h_ws  = 5732.378387768e-6;  
h_wsc = 8385.005556375e-6; 
h_wb  = 13349.334646671e-6;

% Feed water temperature=224C;
T_fw = 224; 

% Tube length
L = 32; %total bundle length
L_b = 21.5;  %this value is fixed to determine accurate heat transfer coefficient in this region.
L_sc = 4.5; %this value is malipulated given that hwsc is known.
L_s = L - L_b - L_sc;

%%%%%Heating circumference%%%%%
R_i = Ri_thermal;
R_o = Ro_thermal;
U_wb  = 2*pi*R_i;
U_ws  = 2*pi*R_i;
U_wsc = 2*pi*R_i;

P_table = 5.0:0.1:6.0; %Pressure;
Ts_avg = (280+318)/2 + 273;
T_table = [];
Hfg_table = [];
hs_table = [];
for PPP = 5.0:0.1:6.0
[T_sat,hf,hg,kf,kg,muf,mug,Prf,Prg] = hsat(PPP);
[dum, hss, dum, dum] = hsh(PPP, Ts_avg);
T_table = [T_table, T_sat];
hs_table = [hs_table, hss];
Hfg_table = [Hfg_table, hg-hf];
end
T_table = (T_table - 273); % Saturated temperature;
Hfg_table = Hfg_table*1e-6;
hs_table = hs_table*1e-6;
a = polyfit(P_table, T_table, 1);
X_5 = a(2); K_5 = a(1);
b = polyfit(P_table, Hfg_table, 1);
X_4 = b(2); K_4 = b(1);
c = polyfit(P_table, hs_table, 1);
dHs_dPs = c(1);

N0 = 8;
% Flow rate on the primary side=560.46kg/s;
W_p = 4717/N/N0;
W_p0 = W_p;
W_p1 = W_p;
W_p2 = W_p;
W_p3 = W_p;
W_p4 = W_p;
W_p5 = W_p;

% Flow rate on the secondary side=62.85kg/s;
W_sec = 502.8/N0/N;
W_s = W_sec;
W_b = W_sec;
W_fw = W_sec;
W_sg = W_sec;
W_s0 = W_sec;
R_i = Ri_thermal;                            
R_o = Ro_thermal; 

%Preparing data matrix;
a_1 = A_p*Cp_p*rho_p/2;   %primary side;
a_2 = A_w*Cp_w*rho_w/2;   %metal;
a_3 = A_s*Cp_s*rho_s/2;   %secondary side;
a_4 = h_pw*pi*R_o/a_1;

a_51 = Cp_p*W_p0/a_1;
a_52 = Cp_p*W_p1/a_1;
a_53 = Cp_p*W_p2/a_1;
a_54 = Cp_p*W_p3/a_1;
a_55 = Cp_p*W_p4/a_1;
a_56 = Cp_p*W_p5/a_1;

a_6 = h_pw*pi*R_o/a_2;
a_7 = h_wb*pi*R_i/a_2;
a_8 = h_ws*pi*R_i/a_2;
a_9 = h_ws*pi*R_i/a_3;

a_11 = A_s*Cp_sc*rho_sc/2; % 
a_12 = h_wsc*pi*R_i/a_2; % 
a_13 = h_wsc*pi*R_i/a_11; % 
a_14 = (Cp_p*rho_p);

% Saturated temperature for 7 Mpa
%c3 = 1000/6.895; % 1Mpa=1000/6.895 lb/inch2
P_s = 5.8;
deltaP = 0.2;
P_sat = (5.8 + deltaP);
%X5=402.94; K5=0.14;  %Tsat~Psat
T_sat = X_5 + K_5*P_sat;
%Tsat=546.6;  %Exit temperature=317C and Degree of superheat is 43.4;
H_fg = X_4 + K_4*P_sat;
sim('HCSG_ss');

T_p1 = Tp10(end);
T_p2 = Tp20(end);
T_p3 = Tp30(end);
T_p4 = Tp40(end);
T_p5 = Tp50(end);
T_p6 = Tp60(end);

T_w1 = Tw10(end);
T_w2 = Tw20(end);
T_w3 = Tw30(end);
T_w4 = Tw40(end);
T_w5 = Tw50(end);
T_w6 = Tw60(end);

T_s1 = Ts10(end);
T_s2 = Ts20(end);
T_sc2 = Tsc20(end);

P_s = P_s;
P_sc = P_s;
%P_sat = P_s;

a_99 = h_ws*pi*R_i/Cp_s;  
a_88 = 2*W_fw/A_s/rho_sc;
a_77 = 2/A_s/rho_sc;
K_sc = -17.615718797133468;       
K_b  = 12.834854076292217;
K_1 = K_5;

P_set = P_s;
Z_ss = 0.76634; % 570K, 60atm; steam expansion coefficient
R = 8.314462E-6; % [MJ/mol-°C]
M_stm       = 0.018;       % [kg/mol]        18.0000       ; % Molar weight of steam          [lbm/lb-mol]

%%%%%%%%%% Pressure Controller Variables %%%%%%%%%%
ts          = 1          ; %
Cst         = 10         ; %
kc          = 5          ; %
ti          = 20         ; %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% BALANCE OF PLANT PARAMETERS %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H_nc        = 2.9193;      % [MJ/kg]         1255.9         ; % nozzle chest enthalpy          [BTU/lbm]
rho_nc      = 27.4877;     % [kg/m^3]           1.716       ; % nozzle chest density           [lbm/ft^3]
w_hp_out    = 158.304;     % [kg/s]      349;%199.0832      ; % mass flow rate out of HP turbin[lbm/s]
rho_rh      = 5.3374;      % [kg/m^3]           0.3332      ; % reheater density               [lbm/ft^3]
H_rh        = 2.8591;      % [MJ/kg]         1230           ; % reheater enthalpy              [BTU/lbm]
w_rh_out    = 33.7926;     % [kg/s]      74.5;%45.3285      ; % output reheat steam flow rate  [lbm/s]
Q_rh        = 5.721113;    % [MJ/s]          5426.2         ; % heat transfer rate in reheater [BTU/s]
w_lp_out    = 99.473;      % [kg/s]    219.3;%125.4166      ; % mass flow rate out of LP turbin[lbm/s]
H_fwh1      = 0.6727;      % [MJ/kg]          289.4181      ; % enthalpy of heater 1           [BTU/lbm]
H_fwh2      = 0.8997;      % [MJ/kg]          387.0721      ; % enthalpy of heater 2           [BTU/lbm]
w_fwh       = 97.0688;     % [kg/s       214;%122.9059      ; % mass flow rate from heater 2 to[lbm/s]
V_nc        = 5.6634;      % [m^3]            200           ; % volume of nozzle chest         [ft^3]
V_rh        = 56.634;      % [m^3]           2000           ; % volume of reheater             [ft^3]
J           = 1e-6;        % [MJ/J]             0.001285067 ; % BTU/ft-lbf
nu_f        = 0.0011;      % [m^3/kg]           0.0184      ; % specific vol. @ 200 psi
A_k2        = 0.2817;      % [m^2]              3.0326      ; % backcalculated from w1=w2
A_k3        = 0.16723;     % [m^2]              1.8         ; % backcalculated from W3 = W2'
h_fwh1      = 0.5927;      % [MJ/kg]     255;%215           ; % heat transfer for heater 1     [BTU/lbm]
h_fwh2      = 0.5927;      % [MJ/kg]     255;%225           ; % heat transfer for heater 2     [BTU/lbm]
tau_fwh1    =         100; % time constant for heater 1
tau_fwh2    =          65; % time constant for heater 2
tau_rh1     =           3; % time constant for reheater flow
tau_rhh     =          10; % time constant for reheater heat transfer
tau_hp      =           2; % time constant for HP turbine
tau_lp      =          10; % time constant for LP turbine
tau_rh2     =          10; % time constant for heater 2 flow
Cp_hr       = 0.0021520;   % [MJ/kg-°C]          0.51399637     ; % Cp of steam @ 531 F and 183 psi. 0.238845896627 BTU/lbm-F = 1 KJ/kg-K
K_b_hp      =      0.1634; % fraction of bleed steam from HP turbine
K_b_lp      =      0.2174; % fraction of bleed steam from LP turbine
eta_hp_isen =       0.861; % isentropic efficiency of HP turbine
eta_lp_isen =       0.861; % isentropic efficiency of LP turbine
H_lp_isen   = 958.4*1054.35e-6/0.45359237; % [MJ/kg]
%R1          =           1/0.01165       ; % ideal gas constant

% Populate look up tables
run('Xsteam/SteamEnthalpy.m');
run('Xsteam/NozzleChestPressure.m');
run('Xsteam/ReheaterPressure.m');

% HP turbine steam enthalpy at isentropic endpoints
term1 = 1080.3*1054.35e-6/0.45359237; % [MJ/kg]
p1    = 200*6894.76;                  % [Pa]
p2    = 1000*6894.76;                 % [Pa]
ce1   = 1.2471e-07;                   % [m^3/kg]
ce2   = 5.3816e-14;                   % [m^6/kg^2]
ce3   = 3.3721e-08;                   % [m^3/kg]

% Condenser Parameters

Pco         =          1735.9         ; % condenser pressure             [psf]
Vco         =          4000           ; % Volume of condenser            [ft^3]
tau         =          10             ; % time constant for condenser sen[s]
Tcwo        =         157.6209       ; % temperature of cooling water ou[F]
Cpco        =         1              ; % Specific heat capacity at const[BTU/lbm-F]
Mcw         =          1000000        ; % cooling water hold up          [lbm]
Wcw0        =         4000;%1397.7         ; % cooling water mass flow rate   [lbm/s]
Tcwi        =         70             ; % temperature of inlet cooling wa[F]
Kco         =          2              ; % condenser PID controller gain
TIco        =         15             ; % PID controller integral time constant
TDco        =         0              ; % PID derivative constant

dLsc=0;
dLb=0;
dLs=0;