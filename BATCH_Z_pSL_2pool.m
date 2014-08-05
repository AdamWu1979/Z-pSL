%   Author: Moritz Zaiss  - m.zaiss@dkfz.de
%   Date: 2014/07/25 
%   Version for cest sources
  
clear all
%% SETUP
 
%pool system parameters
%water pool A
P.R1A=1/3;          % [s^-1] 
P.R2A=2;            % [s^-1]
P.dwA=0; %deltaW_A in [ppm] 

%pool B
P.fB=0.001;         % proton fraction
P.kBA=200;          % [s^-1]
P.dwB=1.9;          % deltaW_B in [ppm} (chemical shift)
P.R2B=30;           % [s^-1]
P.R1B=1;            % [s^-1]
% sequence parameters
P.Zi=1;             % Z initial, in units of thermal M0, Hyperpol.: 10^4                  
P.FREQ=300;         % [MHz]  I use ppm and µT, therefore gamma=267.5153;
P.B1=1;             % [µT]

P.tp=0.1;           % pulse duration [s]
P.n=100;            % pulse number 
P.DC=0.2;           % duty cycle :  saturation time = n * P.tp/DC

P.xZspec= [-5:0.05:5]; % deltaomega [ppm]

Pstart=P;

figure(1), plot(P.xZspec,Z_plsd_2pool(P),'.-') ;   hold on;
set(gca,'XDir','reverse'); xlabel('\Delta\omega [ppm]'); ylabel('Z(\Delta\omega)');
text(min(P.xZspec)/3,0.5,evalc('P'),'FontSize',8)

%% vary a parameter
vary=[ 0.5 1 2 4]; % define value range for variation

for ii=1:numel(vary)
    P=Pstart; % reset previous changes
P.B1=vary(ii); % define which parameter you want to vary
    
figure(2), plot(P.xZspec,Z_plsd_2pool(P),'.-','Color',cl(ii,numel(vary))) ;   hold on;

Pref=P;    Pref.fB=0; %Pref.xZspec=-P.xZspec; %
plot(P.xZspec,Z_plsd_2pool(Pref)-Z_plsd_2pool(P),'.-','Color',cl(ii,numel(vary))) ;   hold on;
% plot(P.xZspec,(P.R1A+P.fC*P.R1C)*(1./Z_cw_2pool(P)-1./Z_cw_2pool(Pref)),'g.-') ;   hold on;

end;
set(gca,'XDir','reverse'); xlabel('\Delta\omega [ppm]'); ylabel('Z(\Delta\omega)');