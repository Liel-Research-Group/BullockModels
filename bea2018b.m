% Bullock et al. (2018b) Probabilistic Models for Tilt
%% INPUTS
% type : flag for model type to use
%   type = 1 for the empirical model for residual tilt
%   type = 2 for the semiempirical model for residual tilt
%   type = 3 for the semiempirical model for peak transient tilt
% S : average foundation settlement in mm
%   only needed for type = 1
% sigS : logarithmic standard deviation of settlement
%   only needed for type = 1
%   if treating settlement as known, leave sigS empty (sigS = [];)
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
% NNS1B : number of non-susceptible layers in the top 1.0B depth of the
%   profile below the foundation
% NS1B : number of susceptible layers in the top 1.0B depth of the profile 
%   below the foundation
% maxHS1B : maximum continuous thickness of susceptible material in the top
%   1.0B of the profile below the foundation
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
% Vgi : peak incremental ground velocity in cm/s
%   the procedure is based on using the median predicted Vgi from the 
%   Bullock et al. (2017) GMPEs; a script for those GMPEs is in the same
%   repository as this script
%% OUTPUTS
% theta : median estimate of residual or peak transient tilt in degrees
% sig : logarithmic standard deviation of tilt
%% CITATION
% Bullock, Z., Dashti, S., Karimi, Z., Liel, A., Porter, K., & Franke, K. 
%    (2018). Probabilistic Models for Residual and Peak Transient Tilt of 
%    Mat-Founded Structures on Liquefiable Soils. Journal of Geotechnical 
%    and Geoenvironmental Engineering, 145(2), 04018108.
%% CODE
function [theta,sig] = bea2018b(type,S,sigS,N160,qc1N,HS,DS,FLPC,NNS1B,NS1B,maxHS1B,q,B,L,Df,h,mass,CAV,Vgi)

DST = DS(1) - HS(1)/2;

DS = DS - Df;

if type == 1
    
    a1 = 0.509;
    a2 = -0.936;
    a3 = -0.102;
    
    sigT = 0.287;
    rho_theta_S = 0.194;
    
    if isempty(sigS) == 0
        sig = sqrt(sigT^2 + a1^2*sigS^2 + 2*a1*rho_theta_S*sigS*sigT);
    else
        sig = sigT;
    end
    
    lnTheta = a1*log(S) + a2*log(B) + a3*DST;
    
    theta = exp(lnTheta);
    
else
   
    alp0 = -4.353;
    alp1 = -0.329;
    alp2 = -0.252;
    alp3 = -0.036;
    alp4 = -0.430;
    alp5 = -0.121;
    alp6 = 0.003;
    alp7 = 0.026;
    alp8 = -0.082;
    alp9 = 0.314;
    alp10 = 0.472;
    alp11 = -0.020;
    alp12 = 0.234;
    alp13 = 0.404;
    
    gam0 = 0.066;
    gam1 = 0.165;
    
    kap0 = 2.383;
    kap1 = 1.491;
    kap2 = -0.168;
    kap3 = -0.327;
    kap4 = 0.087;
    
    sig = 0.548;
    
    if isempty(N160) == 0
        if N160(1) < 17.2
            FLoose = 1;
        else
            FLoose = 0;
        end
    elseif isempty(qc1N) == 0
        if qc1N(1) < 140.2
            FLoose = 1;
        else
            FLoose = 0;
        end
    end
    
    HS1B = 0;
    for i = 1 : length(HS)
        if DS(i) + HS(i)/2 < B
            HS1B = HS1B + HS(i);
        elseif DS(i) - HS(i)/2 < B
            HS1B = HS1B + HS(i) - (DS(i) + HS(i)/2 - B);
        end
    end
    
    lnTheta = alp0 + alp1*log(q) + alp2*log(B)^2 + alp3*L/B + alp4*log(L/B) ...
        + alp5*log(Df) + (alp6*log(Vgi) + alp7*log(CAV))*HS1B + alp8*DS(1) ...
        + alp9*FLoose + alp10*min(mass/10^6,1) + alp11*h/B*mass/10^6 ...
        + alp12*log(Vgi) + alp13*log(CAV);
    
    lnTheta = lnTheta + gam0 + gam1*log(h);
    
    lnTheta = lnTheta + kap0 + kap1*FLPC + kap2*DST + kap3*maxHS1B ...
        + kap4*(NNS1B/NS1B);
    
    theta = exp(lnTheta);
    
    if type == 3
        
        del0 = 2.639;
        del1 = -0.577;
        del2 = -1.285;
        del3 = 0.088;
        
        sigRT = 0.442;
        
        rho_theta_RT = -0.182;
        
        lnRT = del0 + del1*lnTheta + del2*log(B) + del3*HS1B;
        
        lnTheta = lnTheta + lnRT;
        
        theta = exp(lnTheta);
        
        sig = sqrt(sigRT^2 + (1 + del1^2)*sig^2 + 2*del1*rho_theta_RT*sig*sigRT); 
        
    end
    
end

end