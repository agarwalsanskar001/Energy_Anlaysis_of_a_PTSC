clear variables;
%Thermal model of parabolic trough solar collector PTSC

% input parameters of the receiver
D2=0.0158; % meter, inner diameter of the absorber
D3=0.0178; %meter, outer diameter of the absorber
D4=0.057; %meter, inner diameter of the glass envelope
D5=0.060; %meter, outer diameter of the glass envelope
L=3.4; %meter, length of the collector
W=1.1; % width of the collector
In=0.48; % intercept factor
TRS_env=0.9; % transmittance of the glass envelope
ABS_abs=0.9; % absorbance of the absorber
ABS_env=0.02; % Absorptance of the glass envelope

Sigma=5.67*10^-8; % Stefan boltzmann constant


% Input variables
FL=1; % gpm, flow rate of the HTF
u5=5; % m/s, wind speed
T6=296; % K, ambient temperature
T7=T6-8; %K, sky temperature
Month=5; % month's number
n=50; % day's number
T1=45.3; % Inlet temperature
Sstand=13.75; % local time (h)


% Solar radiation code

Sstand1=Sstand*60;
B=(n-1)*360/365;
E=229.2*(0.000075+0.001868*cosd(B)-0.032077*sind(B)-0.014615*cosd(2*B)-0.04089*sind(2*B));
M=4*(60-81)+E;
ST=Sstand1+M;
HA=15*(ST-12*60)/60; % Hour angle
SD=23.45*sind(360*(284+n)/365); %Solar declination
SA=asind(sind(29.1735)*sind(SD)+cosd(29.1735)*cosd(SD)*cosd(HA));
%Solar altitude
ZA=90-SA;%Zenith angle
AOI=ZA; % angle of incidence
azimuth = asind(-sind(HA)*cosd(SD)/ZA);

q_ra0 = 752.08;                   %W/m^2
q_ra = q_ra0*1.1;                 % W/m solar irradiation multiplied by collector width

 
% Optical model

K=3*10^-5*AOI^2-0.0072*AOI+1.2257; % Incidence angle modifier
Endloss=1-0.5/3.6*tand(AOI); % End loss
EEF_abs=In*TRS_env*ABS_abs*K*Endloss; %effective optical efficiency

%Thermal Model

for i=0:0.3:4
    
    T5=294+i;                                   %K estimated outer glass cover temperature

    %Radiation heat transfer from the glass envelope to the sky
        
    E5=0.86;                                    % emissivity of the glass envelope
        
    Rho_D5 = 101325/(287.05*(T6));              % kg/m3, density of the air evaluated at T6
    Mu_D5 = 1.458*10^-6*T6^(3/2)/(T6+110.4);    %kg/m.s, viscosity of the air evaluated at T6
        
    Re_D5=Rho_D5*u5*D5/Mu_D5;                   % Reynolds number at ambient conditions
        
    Pr6 = 0.73025;                              % Prandtl number evaluated at T6
    Pr5 = 0.71;                                 % Prandtl number evaluated at T5
        
    K56=1.5207*10^-11*(T6)^3-4.8574*10^-8*(T6^2)+1.0184*10^-4*(T6)-0.00039333;      % thermal conductivity

    NU_D5= 0.26 * (Re_D5.^(0.6)) *(Pr6.^(0.37)) *((Pr6/Pr5).^0.25);                 % Nusselt number
        
    h56 = K56*NU_D5/D5;             % Convection Heat Transfer Coefficient for air at (T5-T6)/2 (W/m^2-K)

    q_56conv = h56*pi*D5*(T5-T6);   % Convection Heat Transfer from Glass to Ambient
    q_57rad = Sigma*D5*pi*E5*(T5^4-T7^4); %Radiation Heat Transfer
        
    q_loss = q_56conv + q_57rad;        % Total Heat Loss from Glass to Ambient  
        
        
    % Radiation heat transfer from the absorber to the glass envelope
        
    E4 = 0.86;                              % glass envelope emissivity
    Sigma = 5.67*10^-8;                     % Stefan Boltzmann constant
    E3 = 0.0005333*(T1+273.15)-0.0856;      % Absorber select coating emissivity
        
    %conduction heat transfer through the glass envelope

    K45=1.04;                               % glass envelope conductance (Pyrex glass)

    T4=(q_loss + ((T5*(2*pi*K45/log(D5/D4)))/1000) - q_ra*ABS_env*In*Endloss*K)/((2*pi*K45/log(D5/D4))/1000); %inner glass cover temperature
        
    % Convection from HTH to the absorber

    Rho1=1001.1-0.0867*T1-0.0035*T1^2;      % kg/m3, density of the fluid
    Mu1=1.684*10^-3-4.264*10^-5*T1+5.062*10^-7*T1^2-2.244*10^-9*T1^3;       %kg/m.s, viscosity of the fluid
    Cp1=(4.214-2.286*10^-3*T1+4.991*10^-5*T1^2-4.519*10^-7*T1^3+1.857*10^-9*T1^4)*1000;     %J/kg.K, heat capacity of the fluid
        
    K1=0.5636+1.946*10^-3*T1-8.151*10^-6*T1^2;      % W/m.K, conductivity of the fluid
    u1=0.000063*FL/(pi/4*(D2.^2));                  %m/s, HTF fluid speed
    Re_D2=(Rho1*u1*D2)/Mu1;                         % Reynolds number of th HTH fluid
    f2=(1.85*log10(Re_D2)-1.64).^-2;                % friction factor
    Pr1=(Mu1*Cp1)/K1;

    Est_T2=T1+2;
    Ro2=1001.1-0.0867*Est_T2-0.0035*Est_T2^2 ;            % kg/m3, density of the fluid evaluated at T2
    MU2=1.684*10^-3-4.264*10^-5*Est_T2+5.062*10^-7*Est_T2^2-2.244*10^-9*Est_T2^3;        %kg/m.s, viscosity of the fluid evaluated at T2
        
    Cp2=(4.214-2.286*10^-3*Est_T2+4.991*10^-5*Est_T2^2-4.519*10^-7*Est_T2^3+1.857*10^-9*Est_T2^4)*1000 ;    %kJ/kg.K, heat capacity of the fluid evaluated at T2
        
    K2=0.5636+1.946*10^-3*Est_T2-8.151*10^-6*Est_T2^2;    % W/m.K, conductivity of the fluid
        
    Pr2=(MU2*Cp2)/K2;                    % Prandtl number evaluated at T2
        
    NU_D2=((f2/8)*(Re_D2-1000)*Pr1)/(1+12.7*sqrt(f2/8)*(Pr1.^(2/3)-1))*(Pr1/Pr2)^0.11;      %Nusselt number
        
    h1=(NU_D2*K1)/D2;           %W/m2.K, HTF convection heat transfer coefficient
        
    % "conduction heat transfer through the absorber”

    Est_T23= T1+2;                      %C, estimated value
    K23=0.0151*Est_T23+14.5837;         % absorber thermal conductance (321 H) at the average absorber,
        
    %convection from the absorber to the glass envelope

    Kstd=0.02551;               %W/m.K, thermal conductance of annulus gas at standard conditions
    b=1.571;                    %interaction coefficient
    mol_diameter=3.53*10^-8;    % cm , molecular diameter of annulus gas
    lambda=88.67*10^-2;         % m, mean free path between collisions of a molecule
        
    h34=Kstd/((D3/2*log(D4/D3))+b*lambda*(D3/D4+1));        % W/m2.K, convection heat transfer coefficient

        
    syms T2 T3 positive
    
    equ1 = h1*D2*pi*(T2-(T1+273.15))-2*pi*K23/log(D3/D2)*(T3-T2) == 0;
    equ2 = q_ra*(EEF_abs) - pi*D3*h34*(T3-T4) -(Sigma*pi*D3*(T3.^4-T4.^4))/(1/E3+(1-E4)*D3/(E4*D4)) -2*pi*K23/log(D3/D2)*(T3-T2) == 0;

    sol = solve([equ1, equ2], [T2,T3]);

    T2 = double(real(sol.T2));          % inner absorber temperature
    T3 = double(real(sol.T3));          % outer absorber temperature
   

    Qloss=(pi*D3*h34*(T3(1)-T4)+(Sigma*pi*D3*((T3(1)).^4-T4.^4))/(1/E3+(1-E4)*D3/(E4*D4))+q_ra*ABS_env*In*Endloss*K); % heat loss
    Q=(abs(Qloss-q_loss));

    if double(Q) < 1
        break
    end
end
    
q_gain=(h1*D2*pi*(T2(1)-T1-273.15))*L;      %W, heat gain
    
Efficiency=q_gain/q_ra0*100/L/W;            % Efficiency of the collector
    
m=Rho1*u1*(pi/4*(D2.^2));                   % mass flow rate
    
Tout= (((q_ra*(EEF_abs)+q_ra*ABS_env*In*Endloss*K-(Sigma*pi*D3*(T3(1).^4-T4.^4))/(1/E3+(1-E4)*D3/(E4*D4))-pi*D3*h34*(T3(1)-T4))*0.24)/(m*Cp1)+T1);
% outlet temperature

T1=Tout;

