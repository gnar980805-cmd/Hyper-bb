%% Demo
%% Example 1
clearvars
wvs = 532; % wavelength nm
cp = 0.02; % non-water cp, true value without effect of acceptance angle
bbr = 0.005; % backscattering ratio for particles
ssa = 0.7; % single scattering albedo for particles
sigma1 = get_hyperbb_sigmacorr_2026(wvs,cp,bbr,ssa);
% version2 input: wvs, cp, ap, bbr
%sigma1 = get_hyperbb_sigmacorr_v2(wvs,cp,(1-ssa).*cp,bbr); 

% include temperature and salinity
env.S = 35; % default 0 PSU
env.T = 25;% default 20°C
sigma2 = get_hyperbb_sigmacorr_2026(wvs,cp,bbr,ssa,env);
%sigma2 = get_hyperbb_sigmacorr_v2(wvs,cp,(1-ssa).*cp,bbr,env);

%% Example 2
clearvars;
wvs = [430 500 700]; % wavelength(one value or a vector)
cp = [ 0.01 0.1 0.4]; % non-water attenuation
bbr = [0.01 0.02 0.05]; % baclscattering ratio
ssa = [0.1 0.2 0.3]; % single scattering albedo
instr.c_type = 'AC-S'; % AC-S, LISST-VSF
% others.acc_ang = 0.9; % alternatively
sigma = get_hyperbb_sigmacorr_2026(wvs,cp,bbr,ssa,[],instr);
% sigma = get_hyperbb_sigmacorr_v2(wvs,cp,(1-ssa).*cp,bbr,[],instr);


%% Example 3
clearvars
cp = [0 10 20 30 40];
env.T = 25;
env.S = 34;
instr.c_type = 'LISST-VSF'; % AC-S, LISST-VSF

[sigma,data] = get_hyperbb_sigmacorr_2026(430,cp,0.01,0.8,env,instr);
% [sigma,data] = get_hyperbb_sigmacorr_v2(430,cp,0.2.*cp,0.01,env,instr);
plot(cp,sigma)
data

%% Example 4
% in case ACS's absorption doesn't do scattering correction
clearvars
cp = [0 10 20 30 40];
instr.acc_ang = 0.9; 
instr.ang_a = 42;
[sigma,data] = get_hyperbb_sigmacorr_2026(430,cp,0.01,0.1,[],instr);
%[sigma,data] = get_hyperbb_sigmacorr_v2(430,cp,0.9.*cp,0.01,[],instr); 
plot(cp,sigma)
data
