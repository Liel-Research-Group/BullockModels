% Bullock et al. (2018a) Probabilistic Model for Settlement
%% INPUTS
% N160 : vector of SPT blow counts of the susceptible layers
%   **leave empty (N160 = [];) if using CPT data
% qc1N : vector of normalized CPT resistances of the susceptible layers
%   **leave empty (qc1N = [];) if using SPT data
% HS : vector of susceptible layer thicknesses in m
% DS : vector of depths from the ground surface to the center of each 
%   susceptible layer
% FLPC : flag for the presence of a low-permeability cap above the top
%   susceptible layer
%   FLPC = 1 if a low-permeability cap is present
%   FLPC = 0 if NO low-permeability cap is present
% q : foundation bearing pressure in kPa
% B : foundation width in m (shorter foundation plan dimension)
% L : foundation length in m (longer foundation plan dimension)
% Df : foundation embedment depth in m
% h : effective structure height in m
%   for 1 story structures : h = true height of structure
%   for 2+ story structures : h = 0.7 times true height of structure
% mass : inertial mass of the structure in kg
%   the procedure is based on estimating mass based on plan dimensions,
%   number of stories, and construction type
%       mass = 3.41 * (Number of stories) * B * L * rho
%   for wood structures : rho = 51 kg/m3
%   for concrete or masonry structures : rho = 455 kg/m3
%   for steel structures : rho = 242 kg/m3
% CAV : cumulative absolute velocity in cm/s
%   the procedure is based on using the median predicted CAV from the 
%   Bullock et al. (2017) GMPEs; a script for those GMPEs is in the same
%   repository as this script
%% OUTPUTS
% S : predicted average foundation settlement in mm
% sig : logarithmic standard deviation of S
%% CITATION
% Bullock, Z., Karimi, Z., Dashti, S., Porter, K., Liel, A. B., & Franke, 
%    K. W. (2018). A Physics-Informed Semi-Empirical Probabilistic Model 
%    for the Settlement of Shallow-Founded Structures on Liquefiable 
%    Ground. Géotechnique, 1-34.
%% CODE
function [S,sig] = bea2018a(N160,qc1N,HS,DS,FLPC,q,B,L,Df,h,mass,CAV)

NL = length(HS);

DS = DS - Df;

a0 = 1;
if isempty(N160) == 0
    a1 = -0.2174;
else
    a1 = -0.0360;
end
b0 = 0.3026;
b1 = -0.0205;
c0 = 1.3558;
c1 = -0.134;
d0 = -1.3446;
d1 = 0.2303;
d2 = 0.4189;
e0 = -0.8727;
e1 = 0.1137;
e2 = -0.0947;
e3 = -0.2148;
f0 = -0.0137;
f1 = 0.0021;
f2 = 0.1703;
s1 = 0.4973;

k0 = -1.5440;
k1 = 0.0250;
k2 = 0.0295;
k3 = -0.0218;
qc = 61;

fso = zeros(NL);
for i = 1 : NL
    if HS(i) >= 1
        n1 = 0;
        n2 = 0;
        n3 = 0;
        if isempty(N160) == 0
            if N160(i) < 12
                n1 = 1;
            elseif N160(i) < 17
                n2 = 1;
            else
                n3 = 1;
            end
            fso(i) = (n1*a0 + n2*(a0 + a1*(N160(i)-12)) + n3*(a0 + 5*a1))*b0*HS(i)*exp(b1*(max(DS(i),2)^2-4));
        elseif isempty(qc1N) == 0
            if qc1N(i) < 112.4
                n1 = 1;
            elseif qc1N(i) < 140.2
                n2 = 1;
            else
                n3 = 1;
            end
            fso(i) = (n1*a0 + n2*(a0 + a1*(qc1N(i)-112.4)) + n3*(a0 + 27.8*a1))*b0*HS(i)*exp(b1*(max(DS(i),2)^2-4));
        end
    end
end

fsc = FLPC*(c0 + c1*log(CAV));

FSo = sum(fso) + fsc;

fq = (d0 + d1*log(min(CAV,1000)))*log(q)*exp(d2*min(0,B-max(DS(1)-HS(1)/2,2)));
fA = (e0 + e1*log(max(CAV,1500)))*log(B)^2 + e2*L/B + e3*Df;
fh = (f0 + f1*log(min(CAV,1000)))*h^2 + f2*min(mass/10^6,1);

FSt = fq + fA + fh;

logS = FSo + FSt + s1*log(CAV);

[~,critL] = max(fso);

logS = logS + max(k0 + k1*min(HS(critL),12)^2 + k2*min(q,qc) + k3*max(q-qc,0),0);

S = exp(logS);
sig = 0.6746;

end