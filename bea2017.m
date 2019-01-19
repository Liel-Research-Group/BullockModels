% Bullock et al. (2017) Ground Motion Prediction Equations
%% INPUTS
% im : flag for which intensity measure is being predicted
%   im = 1 for cumulative absolute velocity
%   im = 2 for cumulative absolute velocity above 5 cm/s2 threshold
%   im = 3 for standardized cumulative absolute velocity
%   im = 4 for Arias intensity
%   im = 5 for peak incremental ground velocity
% M : moment magnitude of the earthquake scenario
% R : distance to rupture of the earthquake scenario in km
% FR : flag for reverse faulting
%   FR = 1 for reverse faulting
%   FR = 0 otherwise
% FN : flag for normal faulting
%   FN = 1 for normal faulting
%   FN = 0 otherwise
% Z : hypocentral depth in km
% crust : flag for the shallow crustal tectonic environment
%   crust = 1 for shallow crustal earthquakes
%   crust = 0 otherwise
% subduc : flag for the subduction tectonic environment
%   subduc = 1 for interface subduction earthquakes
%   subduc = 2 for intraslab subduction earthquakes
%   subduc = 3 for unknown subduction earthquakes
%   subduc = 0 otherwise
% intra : flag for the intraplate tectonic environment
%   intra = 1 for intraplate earthquakes
%   intra = 0 otherwise
%% OUTPUTS
% Y : intensity measure median estimate in units of cm/s
% sig : logarithmic standard deviation of Y
% ep : probability that the intensity measure exceeds zero (only for CAV5
%   and CAVSTD - empty otherwise)
%% CITATION
% Bullock, Z., Dashti, S., Liel, A., Porter, K., Karimi, Z., & Bradley, B. 
%   (2017). Ground?Motion Prediction Equations for Arias Intensity, 
%   Cumulative Absolute Velocity, and Peak Incremental Ground Velocity for 
%   Rock Sites in Different Tectonic Environments. Bulletin of the 
%   Seismological Society of America, 107(5), 2293-2309.
%% CODE
function [Y,sig,ep] = bea2017(im, M, R, Z, FR, FN, crust, subduc, intra)

ep = [];

if crust == 1
    
    if im == 1 % CAV
        a0 = -5.800;
        a1 = 3.593;
        a2 = -0.231;
        b1 = -1.415;
        b2 = 0.138;
        b3 = -0.007;
        f1 = -0.155;
        f2 = -0.343;
        tau = 0.398;
        phi = 0.331;
    elseif im == 2 % CAV5
        a0 = -9397;
        a1 = 5.356;
        a2 = -0.395;
        b1 = -3.372;
        b2 = 0.381;
        b3 = -0.0110;
        f1 = -0.197;
        f2 = -0.231;
        tau = 0.490;
        phi = 0.530;
        
        alp0 = -7.232;
        alp1 = 3.364;
        bet1 = -1.775;
        bet2 = -0.029;
        bet3 = 0;
        zscore = alp0 + alp1*M + bet1*log(R) + bet2*R + bet3*Z;
        ep = exp(zscore)/(1+exp(zscore));
    elseif im == 3 % CAVSTD
        a0 = -10.836;
        a1 = 5.165;
        a2 = -0.335;
        b1 = -1.563;
        b2 = 0.087;
        b3 = -0.0019;
        f1 = -0.129;
        f2 = -0.273;
        tau = 0.529;
        phi = 0.392;
        
        alp0 = -5.780;
        alp1 = 3.194;
        bet1 = -3.143;
        bet2 = -0.011;
        bet3 = 0;
        zscore = alp0 + alp1*M + bet1*log(R) + bet2*R + bet3*Z;
        ep = exp(zscore)/(1+exp(zscore));
    elseif im == 4 % AI
        a0 = -18.784;
        a1 = 7.758;
        a2 = -0.576;
        b1 = -4.013;
        b2 = 0.435;
        b3 = -0.021;
        f1 = 0.005;
        f2 = -0.580;
        tau = 0.786;
        phi = 0.747;
    elseif im == 5 % Vgi
        a0 = -8.283;
        a1 = 4.480;
        a2 = -0.353;
        b1 = -2.724;
        b2 = 0.294;
        b3 = -0.0064;
        f1 = 0.071;
        f2 = -0.390;
        tau = 0.545;
        phi = 0.452;
    end
    
    lnY = a0 + a1*M + a2*M^2 + (b1 + b2*M)*log(R) + b3*R + f1*FR + f2*FN;
    
    Y = exp(lnY);
    
    sig = sqrt(tau^2+phi^2);
    
elseif intra == 1
    
    if im == 1 % CAV
        a0 = -13.063;
        a1 = 5.078;
        a2 = -0.273;
        b1 = 0.439;
        b2 = -0.145;
        b3 = -0.0047;
        tau = 0.262;
        phi = 0.411;
    elseif im == 2 % CAV5
        a0 = -28.527;
        a1 = 8.034;
        a2 = -0.157;
        b1 = 2.913;
        b2 = -0.825;
        b3 = -0.0089;
        tau = 0.463;
        phi = 0.553;
        
        alp0 = -13.484;
        alp1 = 3.663;
        bet1 = -0.223;
        bet2 = -0.025;
        bet3 = 0;
        zscore = alp0 + alp1*M + bet1*log(R) + bet2*R + bet3*Z;
        ep = exp(zscore)/(1+exp(zscore));        
    elseif im == 3 % CAVSTD
        fprintf('Bullock et al. (2017) does not provide models for standardized CAV in the intraplate tectonic environment.\n');
    elseif im == 4 % AI
        a0 = -33.761;
        a1 = 11.016;
        a2 = -0.717;
        b1 = -0.421;
        b2 = -0.125;
        b3 = -0.0089;
        tau = 0.534;
        phi = 0.325;
    elseif im == 5 % Vgi
        a0 = -3.029;
        a1 = 0.931;
        a2 = -0.040;
        b1 = -0.828;
        b2 = 0.048;
        b3 = -0.0034;
        tau = 0.340;
        phi = 0.483;
    end
    
    lnY = a0 + a1*M + a2*M^2 + (b1 + b2*M)*log(R) + b3*R;
    
    Y = exp(lnY);
    
    sig = sqrt(tau^2+phi^2);
    
elseif subduc ~= 0
    
    if subduc == 2
        delta = 0.00724*10.^(0.507*min(M,8));
    else
        delta = 0.00724*10.^(0.507*min(M,8));
    end
    
    R = sqrt(R^2 + delta^2);
    
    if subduc == 3
        if im == 1 % CAV
            a0 = -4.865;
            a1 = 2.108;
            a2 = 0.036;
            b1 = -0.009;
            b2 = -0.220;
            b3 = 0.0014;
            b4 = 0.0074;
            tau = 0.312;
            phi = 0.302;
        elseif im == 2 % CAV5
            a0 = -0.902;
            a1 = 2.471;
            a2 = 0.131;
            b1 = -2.611;
            b2 = -0.289;
            b3 = 0.0087;
            b4 = 0.0268;
            tau = 0.876;
            phi = 0.869;
        elseif im == 3 % CAVSTD
            a0 = -9.658;
            a1 = 0.705;
            a2 = 0.314;
            b1 = 3.364;
            b2 = -0.702;
            b3 = -0.0035;
            b4 = 0.0191;
            tau = 0.423;
            phi = 0.247;
        elseif im == 4 % AI
            a0 = -15.969;
            a1 = 4.203;
            a2 = 0.092;
            b1 = -0.383;
            b2 = -0.538;
            b3 = 0.0051;
            b4 = 0.0195;
            tau = 0.659;
            phi = 0.683;
        elseif im == 5 % Vgi
            a0 = 1.126;
            a1 = 1.343;
            a2 = -0.062;
            b1 = -2.737;
            b2 = 0.182;
            b3 = 0.0013;
            b4 = 0.0110;
            tau = 0.384;
            phi = 0.421;
        end
    else
        if im == 1 % CAV
            a0 = -3.674;
            a1 = 1.740;
            a2 = 0.098;
            if subduc == 2
                b1 = 0;
                b2 = -0.250;
                b3 = 0;
                b4 = 0.0066;
            elseif subduc == 1
                b1 = 0.381;
                b2 = -0.334;
                b3 = 0.0047;
                b4 = 0;
            end
            tau = 0.319;
            phi = 0.298;
        elseif im == 2 % CAV5
            a0 = -5.796;
            a1 = 1.876;
            a2 = 0.311;
            if subduc == 2
                b1 = 0;
                b2 = -0.572;
                b3 = 0;
                b4 = 0.0161;
            elseif subduc == 1
                b1 = 0.484;
                b2 = -0.699;
                b3 = 0.0067;
                b4 = 0;
            end
            tau = 0.752;
            phi = 0.567;
        elseif im == 3 % CAVSTD
            a0 = 3.684;
            a1 = -0.575;
            a2 = 0.335;
            if subduc == 2
                b1 = 0;
                b2 = -0.357;
                b3 = 0;
                b4 = 0.0089;
            elseif subduc == 1
                b1 = 0.565;
                b2 = -0.503;
                b3 = 0.0071;
                b4 = 0;
            end
            tau = 0.406;
            phi = 0.222;
        elseif im == 4 % AI
            a0 = -15.390;
            a1 = 3.704;
            a2 = 0.201;
            if subduc == 2
                b1 = 0;
                b2 = -0.615;
                b3 = 0;
                b4 = 0.0180;
            elseif subduc == 1
                b1 = 0.692;
                b2 = -0.774;
                b3 = 0.0107;
                b4 = 0;
            end
            tau = 0.582;
            phi = 0.675;
        elseif im == 5 % Vgi
            a0 = -8.735;
            a1 = 2.454;
            a2 = 0.052;
            if subduc == 2
                b1 = 0;
                b2 = -0.289;
                b3 = 0;
                b4 = 0.0095;
            elseif subduc == 1
                b1 = 0.231;
                b2 = -0.354;
                b3 = 0.0060;
                b4 = 0;
            end
            tau = 0.343;
            phi = 0.425;
        end
    end
    
    if im == 2 % CAV5
    
        alp0 = 0.4919;
        alp1 = 3.0830;
        bet1 = -4.0998;
        bet2 = -0.0044;
        bet3 = 0.0278;
        zscore = alp0 + alp1*M + bet1*log(R) + bet2*R + bet3*Z;
        ep = exp(zscore)/(1+exp(zscore));    
    
    elseif im == 3 % CAVSTD
        
        alp0 = 2.3897;
        alp1 = 3.3157;
        bet1 = -5.9813;
        bet2 = 0.0060;
        bet3 = 0.0401;
        zscore = alp0 + alp1*M + bet1*log(R) + bet2*R + bet3*Z;
        ep = exp(zscore)/(1+exp(zscore));    
        
    end
    
    lnY = a0 + a1*M + a2*M^2 + b3*R + b1*log(R) + b2*M*log(R) + b4*Z;
    
    Y = exp(lnY);
    
    sig = sqrt(tau^2 + phi^2);
    
end
