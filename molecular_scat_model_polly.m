%   MOLECULAR_SCAT   Predict molecular scattering
%
%   USAGE:   [mol_return]=molecular_scat(range,lambda,elev,model)
%   
%   range  - in km
%   lambda - lidar wavelength in nm
%   elev   - elevation 90 degree is pointing directly vertically
%          - mol_return is normalised to value at h(0) (default)
%   model  - Model which is used to get the temperature and pressure profile for the calculation of the molecular optical properties 
%  
%

%   Essentially, molecular return C(h) depends on molecular
%   backscatter, B(h), and atmospheric transmittance, T(h),
%         
%      C(h) h^2 = KO(h)  T^2(h) B(h)
%
%   assuming lidar constant K=1 and overlap O(h)=1 (true above 500m).
%   To calculate B(h), need differential Rayleigh scattering cross section
%    
%      B(h) = n(h) [d sigma/d omega] 
%
%   at 180 degrees and atmospheric number density n(h). This varies
%   with temperature and pressure variations and is modelled.
%
%   d sigma/d omega = 5.45 [ lambda./550 ]^-4.09 x 1e-28 cm2sr-1
%
%% NEW OPTICAL FORMULAS REPLACED BY EARLINET STANDARD--> according to Buchholtz 1995
%   To calculate transmittance need atmospheric attenuation coeff.
%
%      T(h) = exp[-2 int{a(h) dh}]
%
%      B(h) = 1.5/4pi a(h) 
%
%  
%    - defined for zenith (ground based) / nadir(space based) 
%        set range accordingly
%    - angle from zenith/nadir in radians

% 

function [mol_return, mol_T, mol_beta, mol_lr]=molecular_scat_model_polly(range, lambda, elev, model)

if nargin~=4
  help molecular_scat;
  return;
end

if isempty(lambda)
  lambda=355;
end

if isempty(elev) 
  elev=pi./2; %nadir/zenith pointing - max angle 90
elseif elev>(pi)
  % assume elev is in degrees, we want radians
  elev = elev./180 .*pi;
end

%!!!
%check whether range or height should be used in interpolation routine...suggest it should be height....

% expects range in m !
if isfield(model,'model_height')
  n = interp1(model.model_height,mean(model.pressure./ model.temperature),range,'linear','extrap');   %287?????-->288.16
else
  % model height is a matrix
  n = zeros(length(model.time),length(range)) * NaN;
  for ii = 1:length(model.time)
    n(ii,:) = interp1(model.height(ii,:),model.pressure(ii,:)./ model.temperature(ii,:),range,'linear','extrap')';
  end
  n = nanmean2(n)';
end
% calculate dz vector
dz=abs(diff(range));
dz=[dz(:);dz(end)]./sin(elev);

%n=interp1(model.height(12,:)./1e3,model.pressure(12,:)./model.temperature(12,:)./287,range,'linear','extrap');

% dsdw=(5.45.*((lambda./550).^(-4.09))).*1e-32; %diff rayleigh xscatter m2 sr-1  hb: where is this eqautin from?
% mol_mass=4.81e-26;               %average molecular mass of air = 4.81e-26 kg
% mol_beta=n./mol_mass.*dsdw;      %molecular backscatter m-1 sr-1
% mol_abs=mol_beta.*4.*pi./1.5;    %molecular absorption m-1
% mol_T=exp(-2.*cumsum(mol_abs.*dz)); %molecular transmittance T=exp(-2 int[mol_abs])
% mol_return=mol_T.*mol_beta; %molecular signal  T^2(h) B(h)


% EARLINET Formulas mostly according to Buchholtz1995
if (lambda==355 || lambda==387)
  roh=0.03010;
end
if (lambda==532 || lambda==607)
  roh=0.02842;
end
if lambda==1064
  roh=0.02730;
end

 nr = 1 + 10 .^(-8) .* (5791817 ./ (238.0185 - (1e3 ./ lambda) .^ 2) + 167909 ./ (57.362 - (1e3 ./ lambda) .^2) ); %refractive index for air, for formula lambda in mum needed, therefore 1e3 
 Ns=2.54743e19; %Number concentration per cm^-3 for standard  conditions
 dsdw=24 .* pi .^ 3 ./ (Ns .^2 .* (lambda .* 10 .^ (-7)) .^ 4) .*((nr .^2 - 1) ./ (nr .^2 + 2)) .^ 2 .* (6 + 3 .* roh) ./ (6 - 7 .* roh); %rayleight cross section per cm^2, therfore conversion from lambda in nm to cm
 mol_abs=dsdw .* Ns .* n ./ 101325.*288.16 .* 100; %per meter, but input in cm scale therfore division by100
 mol_lr=8 .* pi ./ 3 .* (1. + roh ./ 2);
 mol_beta=mol_abs ./ (8 .* pi ./ 3 .* (1. + roh ./ 2)); %molecular backscatter via molecular lidar ratio
 mol_T=exp(-1.*cumsum(mol_abs.*dz)); %molecular transmittance T=exp(-2 int[mol_abs]) 2 removed due to Raman application
 mol_return=mol_T.*mol_beta; %molecular signal  T^2(h) B(h)
 %mol_return=mol_beta;
%
%
%





