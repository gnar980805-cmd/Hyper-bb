function ac_cor = get_ACS_TS_Correction(wa,wc,fname)
% Get ACS Temperature & Salinity correction data
% There are two correction datasets
%   1. Sullivan et al. (2006). The correction dataset is saved in
%   'ACS_T&S_Dependencies.txt'. To apply this correction, use
%   ac_cor = get_ACS_TS_Correction(wa,wc,'S')
%
%   2. Rudiger et al. (2014). Which will be used as default. The correction
%   datasets are saved in two files: 'oe222125093m001.csv',
%   'oe222125093m002.csv'. To apply this correction, use
%   ac_cor = get_ACS_TS_Correction(wa,wc); % or
%   ac_cor = get_ACS_TS_Correction(wa,wc,'R'); 
% output format
%   ac_cor is a structure, with the following fields
%       wl_a dadt datc wl_c dcdt dcds

wa = wa(:);
wc = wc(:);

narginchk(2,3)
if nargin == 2
    % fname = 'ACS_T&S_Dependencies.txt';
    fname = 'R';
end
if fname == 'S' % Sullivan et al correction
    res = readmatrix('ACS_T&S_Dependencies.txt','NumHeaderLines',7);
    wl = res(:,1);
    dadt = interp1(wl,res(:,2),wa);
    dcdt = interp1(wl,res(:,2),wc);
    dcds = interp1(wl,res(:,4),wc);
    dads = interp1(wl,res(:,6),wa);
elseif fname == 'R' % Rudiger et al. correction
    res = readmatrix('PsiT_data_for_appendix.csv','NumHeaderLines',5,'Delimiter',[";"]); % T correction
    wl = res(:,1);
    dadt = interp1(wl,res(:,3),wa);
    dadt(isnan(dadt)) = 0;
    dcdt = interp1(wl,res(:,3),wc);
    dcdt(isnan(dcdt)) = 0;

    res = readmatrix('PsiS_data_for_appendix.csv','NumHeaderLines',5,'Delimiter',[";"]); % S correction
    wl = res(:,1);
    dads = interp1(wl,res(:,5),wa);
    dads(isnan(dads)) = 0;
    dcds = interp1(wl,res(:,5),wc);
    dcds(isnan(dcds)) = 0;
else
    error('Unknown T&S correction method');
end
ac_cor.wva = wa;
ac_cor.dadt = dadt;
ac_cor.dads = dads;
ac_cor.wvc = wc;
ac_cor.dcdt = dcdt;
ac_cor.dcds = dcds;
ac_cor.name = fname;
% ac_cor = [wa dadt dads wc dcdt dcds];
end