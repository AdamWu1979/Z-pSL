% 2pool cw-solution for Z-spectra
% P= struct of system and sequence parameters
% xZspec = frequency offsets in ppm
% returns vector Z(xZspec)

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

function Z=Z_pSL(P)

w_ref=2*pi*P.FREQ;
gamma=267.5153;    % for protons
w1 =P.B1*gamma;
da=(P.xZspec-P.dwA)*w_ref;
theta=atan(w1./da);

%Rex, exchange dependent relaxation in the rotating frame
Rex=Rex_Lorentz(da,w1,(P.xZspec-P.dwB)*w_ref,P.fB,P.kBA,P.R2B);

%Reff, R1rho of pure water
Reff=P.R1A*cos(theta).^2 +P.R2A*sin(theta).^2;

%%
R1rho=Reff+Rex; 
Pz=1;    % for SL           Pz=cos(theta);    % for cw   ; 
Pzeff=1; % for SL           Pzeff=cos(theta); % for cw   ; 

%% pulsed formula 
if P.n<2 P.td=0; else P.td=P.tp*(1/P.DC-1); end; % calc td

Psi = (P.fB-Rex./P.kBA);          %S.Rex ist negativ
        p   = exp(-R1rho*P.tp);                     %S.Rho ist negativ
        r1a = P.R1A+P.fB*P.kBA;
        r1b = P.R1B+P.kBA;
        a   = 1/2*(r1b+r1a+sqrt((r1b-r1a)^2+4*P.fB*P.kBA.^2));
        b   = 1/2*(r1b+r1a-sqrt((r1b-r1a)^2+4*P.fB*P.kBA.^2));
        dAA = 1/(a-b)*(-(b-r1a)*exp(-a*P.td)+(a-r1a)*exp(-b*P.td));
        dAB = P.kBA/(a-b)*(-exp(-a*P.td)+exp(-b*P.td));
        Zss = (1-(1-p.^(-1)).*(Pz.*cos(theta).*P.R1A./(R1rho))-dAA-dAB*P.fB)./ (p.^(-1)-dAA-dAB*Psi);
        xi  = dAA./(1-dAB*Psi).*Pz.*Pzeff;

        Z = (P.Zi-Zss).*xi.^(P.n).*exp((-R1rho*P.tp)*P.n) +Zss;
        
        clear Psi p r1a r1b a b dAA dAB Zss xi

end


%HyperCESTLimit %JCP paper;   assumes (R1B<<kBA, Reff<<R2B)
function Rex=Rex_Hyper_full(da,w1,db,fb,kb,r2b)
ka=kb*fb;
Rex=((ka.*kb.*w1.^2.*((-da+db).^2 + (r2b.*(da.^2 + (ka + kb).^2 + kb.*r2b + w1.^2))./kb))./...
            ((ka + kb).*(db.^2.*w1.^2 + ka.*r2b.*w1.^2) + ...
        (ka + kb).*((da.*db - ka.*r2b).^2 + (db.*ka + da.*(kb + r2b)).^2 + ...
        (ka + kb + r2b).^2.*w1.^2) + (ka + kb + r2b).*(da.^2.*w1.^2 + w1.^4)));
end

%LorentzLimit %NBM paper;    assumes (kAB<<kBA, R1B<<kBA, Reff<<R2B)
function Rex=Rex_Lorentz(da,w1,db,fb,kb,r2b)
ka=kb*fb;

REXMAX= ka.*w1^2./(da.^2+w1.^2).*((da-db).^2 +(da.^2+w1.^2).*r2b./kb + r2b.*(kb+r2b));
GAMMA=2*sqrt( (kb+r2b)./kb.*w1.^2 + (kb+r2b).^2);
Rex=REXMAX./((GAMMA./2).^2+db.^2);
end

   