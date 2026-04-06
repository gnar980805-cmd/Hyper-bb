function [sigma,data]= get_hyperbb_sigmacorr_2026(wvs,cp,bbr,ssa,env,instr)
% get_hyperbb_sigmacorr_2026 - Computing sigma correction factor accounting
% for attenuation of light and potential multiple scattering
% Usage:
%   1. sigma = get_hyperbb_sigmacorr_2026(wvs, cp, bbr，ssa, env)
%   2. sigma = get_hyperbb_sigmacorr_2026(wvs, cp, bbr，ssa, env, instr)
%   Input Arguments
%       wvs     - wavelength [nm]
%       cp      - Non-water attenuation coefficient [m-1]
%       ssa     - single  scattering albedo
%       bbr     - backscattering ratio
%       env.S-salinity [PSU]
%       env.T-Temperature [Degree C]
%       instr.c_type - specify instrument (AC-S, LISST). If c-type is
%       defined, cp represents instrument-measured non-water attenuation
%       coefficient that would include scattering contribution within the
%       acceptance angle
%   or
%       instr.acc_ang - acceptance angle [Degree] of attenuation meter
%       instr.ang_a   - if acs absorption is not corrected for scattering contribution, input 42 [Degree]
%   Output
%       sigma   - sigma factor same size as wvs, cp
%       data    - structs that stores the corrected attenuation, bbr and
%           ssr
%   
%   Usage 1, the default case, where Temperature is 20 degree, salinity
%   is 0, and cp values are considered as true non-water attenuation
%   coefficient
% 
%   Usage 2, if instr doesn't have field c_type or acc_ang, consider cp as
%   true non-water attenuation coefficient
% 
%       instr.c_type = 'ACS';
%       instr.c_type = 'LISST';
%   OR alternatively input acceptance angle
%       instr.acc_ang = [0.9]; % ACS
%
%   Required Functions: get_Water_Absorption_Coef, betasw_HZP2019
%   Author: Lingfeng Zhou, Xiaodong Zhang
%       Division of Marine Science
%       School of Ocean Science and Engineering
%       University of Southern Mississippi
%       Lingfeng.zhou@usm.edu, Xiaodong.zhang@usm.edu

%   To cite:
%   L. Zhou et al., "Angular Weighting Function of Hyper-bb Instrument," 
%   Applied Optics, 2026. (in preparation)

wvs = wvs(:);
cp = cp(:);
bbr = bbr(:);
ssa = ssa(:);
% default values
S = 0;
T = 20;
acc_angle = 0; % acceptance angle
corr_bbr = 1;
corr_ssa = 1;

if nargin == 5 || nargin == 6
    if ~isempty(env)
        if isfield(env,'S') % salinity
            S = env.S;
        end
        if ~isfield(env,'T') % temperature
            T = env.T;
        end
    end
end
if nargin == 6
    if isfield(instr,'c_type') % attenuation meter type
        if strcmpi(instr.c_type,'AC-S') || strcmpi(instr.c_type,'ACS')
            acc_angle = 0.9;
        elseif strcmpi(instr.c_type,'LISST-VSF')  || strcmpi(instr.c_type,'LISSTVSF')
            acc_angle = 0.0387;
        else
            error('Unrecoginized attenuation meter')
        end
    end
    if isfield(instr,'acc_ang') % acceptance angle of attenuation meter
        acc_angle = instr.acc_ang;
    end
    ang_a = 180; % degree
    if isfield(instr,'ang_a') % for AC meter, ang_a = 42 degree
        ang_a = instr.ang_a;
    end
    [corr_bbr,corr_ssa]= bbr_sca_corr(bbr,ssa,acc_angle,ang_a);
end

% correct bbr and ssa
bbr = bbr.*corr_bbr;
ssa = ssa.*corr_ssa;

aw = get_Water_Absorption_Coef(wvs);
ac_cor = get_ACS_TS_Correction(wvs,[],'R');% Rottegers et al. (2014) 
aw_cor = aw + S.*ac_cor.dads+ (T-22.5).*ac_cor.dadt; % Mason:23°C±0.5, Pope and Fry: 22°C±1
awi = aw + (25-22.5).*ac_cor.dadt;
[~,bsw]= betasw_HZP2019(wvs,90,T,S,1);
[~,bswi]= betasw_HZP2019(wvs,90,25,0,1);
cw = aw_cor+bsw(:); % water (in situ)
cwi = awi + bswi(:); % water (mu calibration)

rf_hyperbb = get_rf(bbr,12.5);
rf_cmeter = get_rf(bbr,acc_angle); % apply 
c_true = cp./(1-ssa.*rf_cmeter); % convert measured to true c

cm = c_true.*(1-rf_hyperbb.*ssa); % equivalent c of hyper-bb
sigma = get_sigma(cwi,cm+cw); % use total c
data.bbr_true = bbr;
data.ssa_true = ssa;
data.cp_true = c_true; 
data.rf_hyper = rf_hyperbb;
data.rf_cmeter = rf_cmeter;
end




function sigma = get_sigma(cw,c) % simulated Look-up table
cw = cw(:);
c = c(:);
c_simu = [0 1e-3 1e-2 1e-1 1:50]';
wt_simu = 1.0e-06 .* [
    0.1368
    0.1368
    0.1367
    0.1356
    0.1254
    0.1163
    0.1073
    0.0986
    0.0904
    0.0826
    0.0754
    0.0689
    0.0632
    0.0585
    0.0543
    0.0502
    0.0463
    0.0426
    0.0391
    0.0358
    0.0329
    0.0303
    0.0280
    0.0261
    0.0244
    0.0227
    0.0211
    0.0195
    0.0179
    0.0165
    0.0152
    0.0140
    0.0130
    0.0121
    0.0114
    0.0106
    0.0099
    0.0092
    0.0085
    0.0079
    0.0073
    0.0068
    0.0063
    0.0059
    0.0055
    0.0051
    0.0048
    0.0044
    0.0041
    0.0038
    0.0035
    0.0033
    0.0031
    0.0029];

wt0 = interp1(c_simu,wt_simu,cw);
wt = interp1(c_simu,wt_simu,c);
sigma = wt0./wt;
end

function rf = get_rf(bbr,acc_ang) 
% calculate rf (bf/b) with respect to different acceptance angle and bbr (FF function already applied)
bbr = bbr(:);
rf = zeros(size(bbr));
if acc_ang>0
    [n,m] = FF_phase2_bb_nm(bbr);
    for i = 1:length(bbr)
        fun=@(x)FF_phase2(rad2deg(x),n(i),m(i)).*sin(x);
        rf(i) = integral(fun, 0, deg2rad(acc_ang(1)))*2*pi;
    end
    rf = rf(:);
end

end


function [corr_bbr,corr_ssa]= bbr_sca_corr(bbr,ssa,acc_ang,ang_a)
% [corr_bbr,corr_ssa]= bbr_sca_corr(bbr,ssa,acc_ang_c,acc_ang_a)
% calculate correction factor with respect to measured bba and ssa
% Input:
%       bbr    - measured bbr
%       ssa    - measured ssa
%       acc_ang - acceptance angle of c meter
%           ACS:acc_ang = [0.9]; default
%           ACS:acc_ang = [0.9 42]; if absorption didn't conduct scattering
%               correction
%           LISST: acc_ang = [0.08]
bbr = bbr(:);
ssa = ssa(:);
acc_ang_c = acc_ang;
if nargin<4
    ang_a = 180; % default: if a is corrected
end
% calculate bbr from bbr_measure
bbr_simu =  0.001:0.001:0.5; % real bbr
rf_c_simu = ones(size(bbr_simu));
rf_a_simu = ones(size(bbr_simu));
[n,m] = FF_phase2_bb_nm(bbr_simu);
for i = 1:length(bbr_simu)
    fun=@(x)FF_phase2(rad2deg(x),n(i),m(i)).*sin(x);
    rf_c_simu(i) = integral(fun, 0, deg2rad(acc_ang_c))*2*pi;
    rf_a_simu(i) = integral(fun, 0, deg2rad(ang_a))*2*pi;
end
corr_bbr_simu = rf_a_simu-rf_c_simu;
bbr_mea = bbr_simu./corr_bbr_simu;
% plot(bbr_simu,bbr_mea) % testing
% bbr_true = interp1(bbr_mea,bbr_simu,bbr);
corr_bbr = interp1(bbr_mea,corr_bbr_simu,bbr,'makima');
rf_a = interp1(bbr_mea,rf_a_simu,bbr,"makima"); 
rf_c = interp1(bbr_mea,rf_c_simu,bbr,'makima');
corr_ssa = 1./(rf_a-rf_c+rf_c.*ssa);
end


function [n,m] = FF_phase2_bb_nm(bbr)
% calculate n and m with bb ratio
m_step = 1e-6;
bbr = bbr(:);
m0 = 3.01:m_step:4.7;
n0 = 1.01 +0.1542*(m0-3);
bbr0 = FF_phase2_bb(n0,m0);
[~,ind]=min(abs(bbr0-bbr),[],2);
n = n0(ind);
m=m0(ind);
end

function bb = FF_phase2_bb(n, m)
% calculate backscattering ratio for the Fourier Forand phase function 
% for Junge distributed particles of index: n and of slope m.
% 1 < n < 1.35
% 3.5 < m < 5
% Light scattering by particles in waters, Jonasz and Fournier, 2007, P251

v = (3-m)/2;
d90 = 4/3./((n-1).^2)*(sin(45*pi/180))^2;


tmp1 = 1-d90.^(v+1) - 0.5*(1-d90.^v); % now corrected
tmp2 = (1-d90).*(d90.^v);
bb = 1-tmp1./tmp2;
end


function vsf = FF_phase2(ang, n, m)
% calculate Fourier Forand phase function for Junge distributed particles
% of index: n and of slope m.
% 1 < n < 1.35
% 3.5 < m < 5
% Light scattering by particles in waters, Jonasz and Fournier, 2007, P251

v = (3-m)/2;
d180 = 4/3/((n-1)^2);
sin_ang2 = sin(ang*pi/180/2).^2;
d = d180*sin_ang2;
d1 = 1-d;
dv1 = 1-d.^v;

tmp1 = 4*pi*d1.^2.*d.^v;
tmp2 = v*d1-dv1;
tmp3 = (d.*dv1-v*d1)./sin_ang2;
tmp4 = (1-d180.^v)./(16*pi*(d180-1).*d180.^v).*(3*cos(ang*pi/180).^2-1);
vsf = (tmp2+tmp3)./tmp1+tmp4;
end