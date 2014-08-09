%% analytic Z-spectra
%   Date: 2014/08/01 
%   Version for CEST-sources.de
%   Author: Moritz Zaiss  - m.zaiss@dkfz.de
%   CEST sources  Copyright (C) 2014  Moritz Zaiss
%   **********************************
%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or(at your option) any later version.
%    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%    You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
%   **********************************
%
%   --SHORT  DOC--
%   The parameter struct P contains all system and sequence parameters.
%   Z_cw_2pool(P) calculates the Z-spectrum and returns it as a vector.
%
%   --references--:
%               2-pool-cw:      Zaiss et al. JCP, Zaiss and Bachert NBM      
%               3-pool-cw:      Zaiss et al. NBM    
%               2-pool-plsd SL: Roelloffs et al. NBM 
clear all
%% SETUP
clearvars P Pref Pstart
clc
% setup pool system parameters
%water pool 'a'
P.R1A=1/3;              % longitudinal relaxation rate 1/T1 of pool a  [s^-1]
P.R2A=10;                % transversal relaxation rate 1/T2 of pool a  [s^-1]
P.dwA=0;                % chemical shift of the water pool in [ppm]

% CEST pool 'b'
P.fB=0.0018018;         % proton fraction fB=M0B/M0A, e.g. 50mM creatine in water:  4*50/(2*55.5*1000)=0.0018018;
P.kBA=50;               % exchange rate [s^-1]     % corresponds to creatine at ~ 22�C, pH=6.4 in PBS (Goerke et al.)
P.dwB=1.9;              % chemical shift of the CEST pool in [ppm]
P.R2B=30;               % transversal relaxation rate 1/T2 of pool b  [s^-1]
P.R1B=1;                % longitudinal relaxation rate 1/T1 of pool b  [s^-1]

% sequence parameters
P.Zi=1;                 % Z initial, in units of thermal M0, Hyperpol.: 10^4                  
P.FREQ=300;             % [MHz]  I use ppm and �T, therefore gamma=267.5153;
P.B1=1;                 % [�T]

P.tp=0.1;               % pulse duration [s]
P.n=100;                % pulse number 
P.DC=0.5;               % duty cycle :  saturation time = n * P.tp/DC
                        % DC yields interpulse delay td=P.tp*(1/P.DC-1)

P.xZspec= [-5:0.05:5]; % deltaomega [ppm]

Pstart=P;

figure(1), plot(P.xZspec,Z_pSL(P),'.-') ;   hold on;
set(gca,'XDir','reverse'); xlabel('\Delta\omega [ppm]'); ylabel('Z(\Delta\omega)');
text(0,0,evalc('P'),'FontSize',8)

% PLOT CEST:  Z_cw(Pref)-Z_cw(P) with reference struct Pref.
Pref=P;    
Pref.fB=0;                  % reference: Z of system without CEST pool
%Pref.xZspec=-P.xZspec;     % reference: Z of opposite frequency
plot(P.xZspec,Z_pSL(Pref)-Z_pSL(P),'r.-') ;   hold on; legend({'Z','Z_{ref}-Z'})

%% vary a parameter
vary=[ 0.1:0.1:1]; % define value range for variation

for ii=1:numel(vary)
    P=Pstart; % reset previous changes
P.DC=vary(ii); % define which parameter you want to vary
    
figure(2), plot(P.xZspec,Z_pSL(P),'.-','Color',cl(ii,numel(vary))) ;   hold on;

Pref=P;    Pref.fB=0; %Pref.xZspec=-P.xZspec; %
plot(P.xZspec,Z_pSL(Pref)-Z_pSL(P),'.-','Color',cl(ii,numel(vary))) ;   hold on;
% plot(P.xZspec,(P.R1A+P.fC*P.R1C)*(1./Z_cw_2pool(P)-1./Z_cw_2pool(Pref)),'g.-') ;   hold on;

end;
set(gca,'XDir','reverse'); xlabel('\Delta\omega [ppm]'); ylabel('Z(\Delta\omega)');