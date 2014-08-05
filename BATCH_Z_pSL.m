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