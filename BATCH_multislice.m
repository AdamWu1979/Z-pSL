%% analytic Z-spectra
%   Date: 2021/08/01 
%   Version for CEST-sources.de
%   Author: Moritz Zaiss  Jan-Rüdiger Schüre
%   CEST sources - Copyright (C) 2021  Moritz Zaiss
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
%   Thsi code simulates a magnetization history, which occurs for example
%   in multi-splice CEST with interleaved saturation ( Ellingson et al,
%   etc. 
%
%   --references--:
%   2-pool-cw:      Zaiss and Bachert. NBM 2013;26(5):507–18. doi:10.1002/nbm.2887   and Zaiss et al JCP 2012;136:144106. doi:10.1063/1.3701178     
%   3-pool-cw:      Zaiss et al. NBM 2015 Feb;28(2):217-30. doi: 10.1002/nbm.3237.   
%   2-pool-pulsd SL: Roeloffs et al. NBM 2014; 28, 40–53, doi: 10.1002/nbm.3192. 

clear all


%water pool 'a'
P.R1A=1/3;              % longitudinal relaxation rate 1/T1 of pool a  [s^-1]
P.R2A=1/0.8;            % transversal relaxation rate 1/T2 of pool a  [s^-1]
P.dwA=0;                % chemical shift of the water pool in [ppm]

% CEST pool 'b'
P.fB=0.65/1000;         % proton fraction fB=M0B/M0A, e.g. 50mM creatine in water:  4*50/(2*55.5*1000)=0.0018018;
P.kBA=30;               % exchange rate [s^-1]     % corresponds to creatine at ~ 22ï¿½C, pH=6.4 in PBS (Goerke et al.)
P.dwB=3.5;              % chemical shift of the CEST pool in [ppm]
P.R2B=1/0.033;               % transversal relaxation rate 1/T2 of pool b  [s^-1]
P.R1B=1/0.77;                % longitudinal relaxation rate 1/T1 of pool b  [s^-1]

% sequence parameters
P.Zi=1;                 % Z initial, in units of thermal M0, Hyperpol.: 10^4                  
P.FREQ=123.26;            % [MHz]  I use ppm and 3T, therefore gamma=267.5153;
%P.FREQ=300;             % [MHz]  I use ppm and 7T, therefore gamma=267.5153;



P.B1=sqrt(1*0.5)
P.tp=0.25;            % pulse duration [   
P.Trec=0.25;             
P.n=1;                % pulse number 
P.DC=0.5;             % duty cycle :  saturation time = n * P.tp/DC
P.nslices=16;

clear zstick
xZspec= [-8:0.4:8]; % deltaomega [ppm]
help1=round(numel(xZspec)/2)+1;  
mtrrange=xZspec(help1:end);    % X-range MTRasym
Z=1;
for fo=1:numel(xZspec)
   
    P.xZspec= xZspec(fo); % deltaomega [ppm]

        for ii=1:P.nslices
                P.Zi=(Z-1)*exp(-P.R1A*P.Trec)+1;  % Mit Relaxtionsterm  --> gepulste Gl.
                Z=Z_pSL(P);
                zstick(fo,ii)=Z;
             
end
end

cj=jet(P.nslices);
h=figure
for l=1:P.nslices
plot(xZspec,zstick(:,l),'.-','Color',cj(l,:),'Displayname', sprintf('slice %d, = %d sat modules',l,l)); hold on; title('The Z-spectrum moves due to history of the magnetozation')
% plot(mtrrange,flipud(zstick(1:help1-2,l))-zstick(help1:end,l),'.-','Color',cj(l,:)) ; hold on;
end
legend show;
set(gca,'XDir','reverse'); xlabel('\Delta\omega [ppm]'); ylabel('Z(\Delta\omega)');
