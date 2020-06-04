function [] = target_cat_only_lambdaatt_cloud_new(thedate, site, iswrite)

if nargin < 1
  help target_cat_only;
  return;
end

clear plot


%=========set directories===
%Pathes and output file names
cat_path         = [get_data_dir '/' site '/products/classification/' thedate(1:4) '/'];
path_calibrated  = [get_data_dir '/' site '/calibrated/pollyxt/' thedate(1:4) '/'];
path_aclass      = [get_data_dir '/' site '/products/mwl_aclass/' thedate(1:4) '/'];
%cat_path         = [get_data_dir '/experimental/' site '/products/classification/' thedate(1:4) '/'];
%path_calibrated  = [get_data_dir '/experimental/' site '/calibrated/pollyxt/' thedate(1:4) '/'];
%path_aclass      = [get_data_dir '/experimental/' site '/products/mwl_aclass/' thedate(1:4) '/'];

%=======set file names======
file_calibrated  = [path_calibrated  thedate, '_' site '_pollyxt.nc'];
file_aclass      = [path_aclass thedate '_' site '_mwl_aclass.nc'];
not_show_bad_values=0;

%================read data====================
%Load and concatenate the netcdf files
    disp(['Loading Polly data from file : ' file_calibrated]);
	filename        = [file_calibrated];
	[data,attribute] = load_nc_struct(filename);
    min_heightbin = 1;   

if isfield(data,'quality_mask_1064') 
  data.att_beta_0355(data.quality_mask_0355 > 0) = NaN;
  data.att_beta_0532(data.quality_mask_0532 > 0) = NaN;
  data.att_beta_1064(data.quality_mask_1064 > 0) = NaN;
  data.att_beta_0387(data.quality_mask_0387 > 0) = NaN;
  data.att_beta_0407(data.quality_mask_0407 > 0) = NaN;
  data.att_beta_0607(data.quality_mask_0607 > 0) = NaN;
  data.volume_depolarization(data.quality_mask_volume_depolarization > 0) = NaN;

end

% switch lower(site)
%   case 'leipzig'
%     minheight  = 500;
%   case 'krauthausen'
%     minheight  = 500;
%   case 'melpitz'    
%     minheight  = 500;
% end
%     
%====================Preparation for target mask






%======Correct incomplete overlap====================

%overlap_0355_ascii = dlmread('pollyXT_ift_20130418-0000-0145-55sr-355.o34','\t',1,1);
%overlap_0532_ascii = dlmread('pollyXT_ift_20130418-0000-0145-55sr-532.o34','\t',1,1);
%starts reading at row offset R1 and column offset C1. For example, the offsets R1=0, C1=0 specify the first value in the file. 
%To specify row and column offsets without specifying a delimiter, use an empty string as a placeholder, for example, M = dlmread(filename,'',2,1)
%overlap_0355 = 1./(imresize(overlap_0355_ascii(:,1),(size(data.range))));
%overlap_0532 = 1./(imresize(overlap_0532_ascii(:,1),(size(data.range))));
%data.att_beta_0355=data.att_beta_0355.*(overlap_0355*ones(size(data.att_beta_0355),1)')';
%data.att_beta_0387=data.att_beta_0387.*(overlap_0355*ones(size(data.att_beta_0387),1)')';
%data.att_beta_0532=data.att_beta_0532.*(overlap_0532*ones(size(data.att_beta_0532),1)')';
%data.att_beta_0607=data.att_beta_0607.*(overlap_0532*ones(size(data.att_beta_0607),1)')';
%data.att_beta_1064=data.att_beta_1064.*(overlap_0532*ones(size(data.att_beta_1064),1)')';

%======Interpolate to surface==
%attention, when ever signal ratios are implemented this interpolation must
%be rueckgaengig gemacht erden or the respective variables renames so that
%for the signal ratios no interpolation is performed


interpolheight=500;
index = find(data.height(:,1) < interpolheight); 
maxindex=(max(index)+1);
interpol_before=1
if interpol_before==1
    data.att_beta_1064(:,index(:,1))=data.att_beta_1064(:,maxindex)*(ones(1,(maxindex-1)));
    data.att_beta_0532(:,index(:,1))=data.att_beta_0532(:,maxindex)*(ones(1,(maxindex-1)));
    data.att_beta_0355(:,index(:,1))=data.att_beta_0355(:,maxindex)*(ones(1,(maxindex-1)));
    data.att_beta_0387(:,index(:,1))=data.att_beta_0387(:,maxindex)*(ones(1,(maxindex-1)));
    data.att_beta_0607(:,index(:,1))=data.att_beta_0607(:,maxindex)*(ones(1,(maxindex-1)));
end
% old, but loops ar enot effcient in matlab so removed
%by statement above
%for i=1:maxindex 
%    aa=i
%data.att_beta_1064(:,i)=data.att_beta_1064(:,maxindex);
%end



%======smooth temporally====================
smooth_t=10%20 % xx profiles are summed up , i.e. 10 means in case of polly that 10 30 sec profiles are summed so that the temporal resolution is 5 min.
smooth_v=1%9

if smooth_t>1
  data.att_beta_1064   = smooth_temp(data.att_beta_1064,smooth_t)./smooth_t;
  data.att_beta_0532   = smooth_temp(data.att_beta_0532,smooth_t)./smooth_t;
  data.att_beta_0355   = smooth_temp(data.att_beta_0355,smooth_t)./smooth_t;
  data.att_beta_0387   = smooth_temp(data.att_beta_0387,smooth_t)./smooth_t;
  data.att_beta_0607   = smooth_temp(data.att_beta_0607,smooth_t)./smooth_t;
  data.volume_depolarization = smooth_temp(data.volume_depolarization,smooth_t)./smooth_t;
  data.time = smooth_temp_axis(data.time,smooth_t);
end

if smooth_v>1
 
 data.att_beta_1064   = smooth_vert(data.att_beta_1064,smooth_v)./smooth_v;
 data.att_beta_0532   = smooth_vert(data.att_beta_0532,smooth_v)./smooth_v;
 data.att_beta_0355   = smooth_vert(data.att_beta_0355,smooth_v)./smooth_v;
 data.att_beta_0387   = smooth_vert(data.att_beta_0387,smooth_v)./smooth_v;
 data.att_beta_0607   = smooth_vert(data.att_beta_0607,smooth_v)./smooth_v;
 data.volume_depolarization = smooth_vert(data.volume_depolarization,smooth_v)./smooth_v;
 data.range = smooth_vert_axis(data.range,smooth_v);
 data.height = smooth_vert_axis(data.height,smooth_v);
end

%=============================
% Noise reduction  
%=============================
% Can improve this by looking at running mean of data..
% create index array [ntimes,nheights] where all noisy pixels are set to 1 and all valid pixels are set to 0

%--
% despeckle
% despeckle remove quite a lot of valid pixels
%--
%index = despeckle(data.beta_0355);
%data.beta_0355 = remove_inner_pixels(data.beta, index,'despeckle');
%index = despeckle(data.depolarization);
%data.depolarization = remove_inner_pixels(data.depolarization, index,'despeckle');
%=============================
% Finished Noise reduction
%=============================


lidar_ratio =55;%50
ang_exp =1;%1.34;
dz=abs(data.range(2)-data.range(1))
%=============================
% Molecular and particle scattering
%============================
%%---- theoretical molecular scattering profiles ----%%
model = load_model(thedate, site); %load model
%1064
 [mol_return, mol_T, mol_beta] = molecular_scat_model_polly(data.range, 1064, [], model);  % new EARLINET formulas, mol values at 532 nm
 atten_molecular_1064 = ones(length(data.time),1) * (mol_T'); % generate  field of molecular transmission
 molecular_beta_1064 = ones(length(data.time),1) * (mol_beta'); % generate  field of molecular backscattering
 data.corr_beta_1064   = data.att_beta_1064   ./ (atten_molecular_1064.*atten_molecular_1064);  % correct calibrated bsc for molecular attenuation, output m^-1 sr-^1
 data.par_beta_1064  =  data.corr_beta_1064 - molecular_beta_1064;
 index = find(data.par_beta_1064 < 0);
 data.par_beta_1064(index) = 0;
  if interpol_before==0
    data.par_beta_1064(:,index(:,1))=data.par_beta_1064(:,maxindex)*(ones(1,(maxindex-1)));
 end
 %calculate particle beta with neglected particle attenuation
 data.par_ex_1064 = (data.par_beta_1064.*lidar_ratio);
 atten_particle_1064 =  exp(-1.*nancumsum(data.par_ex_1064,2).*dz);
 data.par_beta_1064 =  data.att_beta_1064   ./ (atten_molecular_1064.*atten_molecular_1064.*atten_particle_1064.*atten_particle_1064)-molecular_beta_1064;%new!!!!
 time=(data.time)';
 %size(time);
  
%355
[mol_return, mol_T, mol_beta] = molecular_scat_model_polly(data.range, 355, [], model);  % new EARLINET formulas, mol values at 355 nm
atten_molecular_0355 = ones(length(data.time),1) * (mol_T');
molecular_beta_355 = ones(length(data.time),1) * (mol_beta'); % generate  field of molecular backscattering
data.corr_beta_0355   = data.att_beta_0355   ./ (atten_molecular_0355.*atten_molecular_0355); % correct calibrated bsc for molecular attenuation, output m^-1 sr-^1 
data.par_beta_0355  =    data.corr_beta_0355 - molecular_beta_355;
index = find(data.par_beta_0355 < 0);
data.par_beta_0355(index) = 0;
if interpol_before==0
    data.par_beta_0355(:,index(:,1))=data.par_beta_0355(:,maxindex)*(ones(1,(maxindex-1)));
end
%calculate particle beta with neglected particle attenuation
data.par_ex_0355 = (data.par_beta_0355.*lidar_ratio);
atten_particle_0355 =  exp(-1.*nancumsum(data.par_ex_0355,2).*dz);
%plot(data.par_ex_0355(200,1:150),data.range(1:150))
%data.par_beta_0355 = data.par_beta_0355./(atten_particle_0355.*atten_particle_0355);    %old
data.par_beta_0355 =  data.att_beta_0355   ./ (atten_molecular_0355.*atten_molecular_0355.*atten_particle_0355.*atten_particle_0355)-molecular_beta_355;%new!!!!
%387
[mol_return, mol_T, mol_beta] = molecular_scat_model_polly(data.range, 607, [], model); 
atten_molecular_0387 = ones(length(data.time),1) * (mol_T');  % generate  field of molecular transmission
molecular_beta_0387 = ones(length(data.time),1) * (mol_beta');

%532
 [mol_return, mol_T, mol_beta] = molecular_scat_model_polly(data.range, 532, [], model);  % new EARLINET formulas, mol values at 532 nm
 atten_molecular_0532 = ones(length(data.time),1) * (mol_T');  % generate  field of molecular transmission
%atten_particle_0532 =  exp(-1.*nancumsum(data.par_ex_0532,2).*dz);
 molecular_beta_532 = ones(length(data.time),1) * (mol_beta'); % generate  field of molecular backscattering
 data.corr_beta_0532   = data.att_beta_0532   ./ (atten_molecular_0532.*atten_molecular_0532);  % correct calibrated bsc for molecular attenuation, output m^-1 sr-^1
 data.par_beta_0532  = data.corr_beta_0532 - molecular_beta_532;
 index = find(data.par_beta_0532 < 0);
 data.par_beta_0532(index) = 0;
  if interpol_before==0
    data.par_beta_0532(:,index(:,1))=data.par_beta_0532(:,maxindex)*(ones(1,(maxindex-1)));
  end%calculate particle beta with neglected particle attenuation
 data.par_ex_0532 = (data.par_beta_0532.*lidar_ratio);
 atten_particle_0532 =  exp(-1.*nancumsum(data.par_ex_0532,2).*dz);
% data.par_beta_0532 =
% data.par_beta_0532./(atten_particle_0532.*atten_particle_0532); old
data.par_beta_0532 =  data.att_beta_0532   ./ (atten_molecular_0532.*atten_molecular_0532.*atten_particle_0532.*atten_particle_0532)-molecular_beta_532;%new!!!!
 
 
 
 
 %Particle depolarization
 dmol=0.0072%0.0125%0.0072;%53; %for 1 nm filters
 data.pardepol = (data.volume_depolarization + 1) ./ (molecular_beta_532 .* (dmol*ones(length(data.time), length(data.range)) - data.volume_depolarization) ./ (data.par_beta_0532 .* (1 + dmol) ) +1 )-1; %"particle depol" at 532 nm

%607
 [mol_return, mol_T, mol_beta] = molecular_scat_model_polly(data.range, 607, [], model); 
 atten_molecular_0607 = ones(length(data.time),1) * (mol_T');  % generate  field of molecular transmission
 molecular_beta_0607 = ones(length(data.time),1) * (mol_beta');

 AOD_0355 = (trapz(data.par_ex_0355,2)).*dz;
 AOD_0532 = (trapz(data.par_ex_0532,2)).*dz;
 AOD_1064 = (trapz(data.par_ex_1064,2)).*dz;
 save('aodvrgl.txt','time','AOD_0355','AOD_0532','AOD_1064','-ascii','-double', '-tabs');
%Angström exponents

data.ang355532  =log(data.par_beta_0532 ./ data.par_beta_0355) ./ log(355/532);
%index = find((data.par_beta_0355<1e-7) | (data.par_beta_0532<1e-7));
%data.ang355532(index) = nan;

data.ang5321064 =log(data.par_beta_1064 ./ data.par_beta_0532) ./ log(532/1064);
%index = find((data.par_beta_0532<1e-7) | (data.par_beta_1064<1e-7));
%data.ang5321064(index) = nan;

data.ang3551064 =log(data.par_beta_1064 ./ data.par_beta_0355) ./ log(355/1064);
%index = find((data.par_beta_0355<1e-7) | (data.par_beta_1064<1e-7));
%data.ang3551064(index) = nan;

%=============================
% Finished Molecular and Particle Scattering
%=============================

%======================================
%Mass concentration
%======================================
% 
% depol_d = 0.31;
% depol_nd= 0.05;
% rho_d = 2.6;
% conversion_factor_dust= 0.65E-6;
% S_d = 50;
% data.dust_beta_0532=zeros(size(data.par_beta_0532));
% data.dust_beta_0532  = data.par_beta_0532.*(((data.pardepol-depol_nd).*(1+depol_d))./((depol_d-depol_nd).*(1+data.pardepol)));
% index2=find(data.pardepol <= depol_nd);
% data.dust_beta_0532(index2)=0;
% data.dust_mass_conc=1000000000000.*rho_d.*conversion_factor_dust.*data.dust_beta_0532*S_d;
% fff=size(data.dust_mass_conc)
% jjj=data.dust_mass_conc(1,25)
% 

%=============================================================
%Read Cloudnet classification
%=============================================================

% filedate = sprintf('%04d%02d%02d', str2double(thedate(1:4)), str2double(thedate(5:6)), str2double(thedate(7:8)));
% filename        = [cat_path filedate '_' site '_classification.nc']
% 
% 		[cloudnet cloudnet_att] = load_nc_struct(filename);
%         aaaa=cloudnet.target_classification;
%         aaab=cloudnet.time;
%         aaac=cloudnet.height;
		


%=============================
% Aerosol Typing
%=============================
% Flag data - cloud, aerosol, molecular, noise

% 0: Noise
% 1: molecular
% 2: cloud/aerosol but not typed 
% 3: cloud (liquid, near base before ms kicks in)
% 4: aerosol
% 5: ice containing cloud 
% 6: dust

%start mask with everything is noise
mask = zeros(size(data.att_beta_0355));

% signal present in 355--> molecular
index = find(~isnan(data.att_beta_0355 ));
size(index)
mask(index) = 1;
%-----------------------------------------------------------------------
%new mask
% beta 1064 "sees" something
 index = find(data.par_beta_1064 > 1e-8)   ; % m^-1 sr-^1  %untyped aerosol
 mask(index) =2;
 %size(index)
 % Aerosol
 index = find(data.par_beta_1064 > 2e-7  & data.ang5321064 >= 0.75 & data.pardepol < 0.07);% | isnan(data.pardepol) ));%%small aerosol vorher 1.5 und ang355532
 mask(index) = 3;
 index = find(data.par_beta_1064 > 2e-7 & data.par_beta_0532 > 2e-7  & data.pardepol < 0.2 & data.pardepol >= 0.07);% & abs(data.ang5321064) <=0.5); %polluted dust
 mask(index) = 5;
 index = find(data.par_beta_1064 > 2e-7 & data.par_beta_0532 > 2e-7  & data.pardepol >= 0.2);% & abs(data.ang5321064) <= 0.5 );%& data.pardepol <= 0.35); %dust
 mask(index) = 6;
 index = find(data.par_beta_1064 > 2e-7 &  abs(data.ang5321064) <.75 & data.pardepol < 0.07);% | isnan(data.pardepol)));%spherical,large aerosol
 mask(index) = 4;

% Liquid layers have low depolarization but high backscatter 
%index = find(data.par_beta_1064 > 2e-5); %& abs(data.ang5321064) <= 0.5); %untyped cloud
%mask(index) = 7; 
%index = find(data.par_beta_1064 > 2e-5 & data.pardepol <= 0.05); %most prob drops
%mask(index) = 9; 
%index = find(data.par_beta_1064 > 5e-5 & data.pardepol <= 0.05 & abs(data.ang5321064) <= 0.5); %drops
    droplet_bit = get_droplet_bit(data.height, data.par_beta_1064);
    mask(find(droplet_bit==1)) = 7;
    mask(find(droplet_bit==1 & data.pardepol <= 0.05)) = 9;
    mask(find(droplet_bit==1 & data.pardepol <= 0.05 & abs(data.ang5321064) <= 0.5)) = 8;

%Ice
index = find(data.par_beta_1064 > 2e-7 & data.par_beta_0532 > 2e-7 & data.volume_depolarization >= 0.3); %most prob ice
mask(index) = 11;
index = find(data.par_beta_1064 > 2e-7 & data.par_beta_0532 > 2e-7 & data.pardepol >= 0.45); %& abs(data.ang5321064) <= 0.8); %ice vorher 1e-6 threshold
mask(index) = 10;

%data.flag = mask;

%post-processing
maskx=size(mask,1);
masky=size(mask,2);
%if liquid cloud found, set everything above which ias not cloud anymore to
%0--> no siognal
 for i=1:maskx
   index=find(((mask(i,:)>6)&(mask(i,:)<10)),1,'first');
   index2=find((mask(i,(index+1):masky)<7 | mask(i,(index+1):masky)>9));  %changed from <9 to >9, +1 added
   mask(i,((index+1)+index2-1))=0;%+1 added
 end    
%for i=1:maskx
   %  index=find(((mask(i,:)>6)&(mask(i,:)<10)),1,'first');
  %   index2=find((mask(i,(index+1):masky)<8 | mask(i,(index+1):masky)>9),1,'first');  %changed from <9 to >9, +1 added
 %    mask(i,((index+1)+index2-1):masky)=0;%+1 added
%end 
%set mask to no signal below interpol heigth
 index = find(data.height<interpolheight);
 for i=1:maskx
     mask(i,index) = 0;
     data.par_beta_0355(i,index)  = NaN;
     data.par_beta_0532(i,index)  = NaN;
     data.par_beta_1064(i,index)  = NaN;
  end

%=============================
% NetCDF Output
%(active, if argument 'iswrite'=1)
%=============================
if iswrite
%%reduce data arrays to max_heightbins
max_height    = 10000; %m

%%-- Correct missing values --%%
  missing_value = -999;

  %%--- setup netcdf file ---%%
  if ~isempty(data.time)
    if ~exist(path_aclass,'dir')
      status=mkdir(path_aclass);
     if ~status; disp(['Could not create directory in ' path_aclass]); return; end;
    end
    disp(['Storing mwl target classification to' file_aclass]);

    g = netcdf(file_aclass,'clobber');
    %%--- write dimensions to netcdf file ---%%
    g('time')  = length(data.time);
    g('range') = length(data.range)
    g('channels') = 3; %added

    %%--- write variables and attributes to netcdf file ---%%
    g{'time'} = ncfloat('time');
    g{'time'}(:) = data.time;
    g{'time'}.units = ['hours since ' thedate ' 00:00:00 +00:00'] ;
    g{'time'}.long_name = 'Decimal hours from midnight UTC to the middle of each ray' ;
    g{'time'}.axis = 'T' ;
    g{'time'}.missing_value = missing_value;
    g{'time'}.FillValue_ = missing_value;

    g{'range'} = ncfloat('range');
    g{'range'}(:) = data.range;
    g{'range'}.units = 'm' ;
    g{'range'}.long_name = 'Range from lidar to centre of range bin' ;

    g{'height'} = ncfloat('range');
    g{'height'}(:) = data.range*cos(5.*pi./180);
    g{'height'}.units = 'm' ;
    g{'height'}.long_name = 'Height of centre of range bin above lidar' ;
    g{'height'}.axis = 'Z' ;

    g{'latitude'} = ncfloat;
    g{'latitude'}(:) = data.latitude;
    g{'latitude'}.units = 'degrees_north' ;
    g{'latitude'}.long_name = 'Latitude of lidar' ;
    g{'latitude'}.standard_name = 'latitude' ;

    g{'longitude'} = ncfloat;
    g{'longitude'}(:) = data.longitude;
    g{'longitude'}.units = 'degrees_east' ;
    g{'longitude'}.long_name = 'Longitude of lidar' ;
    g{'longitude'}.standard_name = 'longitude' ;

    g{'altitude'} = ncfloat;
    g{'altitude'}(:) = data.altitude;
    g{'altitude'}.units = 'm' ;
    g{'altitude'}.long_name = 'Height of instrument above mean sea level' ;

    g{'lidar_beam_divergence'} = ncfloat;
    g{'lidar_beam_divergence'}(:) = 0.2;
    g{'lidar_beam_divergence'}.long_name = 'Lidar laser beam divergence';
    g{'lidar_beam_divergence'}.units = 'mrad';

    g{'lidar_telescope_field_of_view'} = ncfloat;
    g{'lidar_telescope_field_of_view'}(:) = 2;
    g{'lidar_telescope_field_of_view'}.long_name = 'Lidar telescope field of view';
    g{'lidar_telescope_field_of_view'}.units = 'mrad';
    
    g{'angstroem'} = ncfloat;
    g{'angstroem'}(:) = ang_exp;
    g{'angstroem'}.units = '1' ;
    g{'angstroem'}.long_name = 'Angstroem exponent used for the calculation of the quasi backscatter coefficients' ;

    g{'lidar_ratio'} = ncfloat;
    g{'lidar_ratio'}(:) = lidar_ratio;
    g{'lidar_ratio'}.units = 'sr' ;
    g{'lidar_ratio'}.long_name = 'lidar ratio at 1064 nm used to obtain the quasi particle exctinction coefficient' ;

    
    g{'flag'} = ncbyte('time','range');
    g{'flag'}(:) = mask;
    g{'flag'}.long_name = 'Lidar target classification';
    g{'flag'}.plot_range = ncbyte([0 11]);
    g{'flag'}.legend_key_red  =[255 , 0.9*255 , 0.6*255 , 221 , 231 , 136 , 0 , 120 , 058 , 180 , 017 , 134]/255;
    g{'flag'}.legend_key_green=[255 , 0.9*255 , 0.6*255 , 204 , 109 , 034 , 0, 028 , 137 , 221 , 119 , 187]/255;
    g{'flag'}.legend_key_blue =[255 , 0.9*255 , 0.6*255 , 119 , 046 , 000 , 0 , 129 , 201 , 247 , 051 , 106]/255;

    g{'flag'}.definition= ...             
     ['0: No signal' 10 ...
      '1: Clean atmosphere' 10 ...
      '2: Non-typed particles/low conc.' 10 ...
      '3: Aerosol: small' 10 ...
      '4: Aerosol: large, spherical' 10 ...
      '5: Aerosol: mixture, partly non-spherical' 10 ...
      '6: Aerosol: large, non-spherical' 10 ...
      '7: Cloud: non-typed' 10 ...
      '8: Cloud: water droplets' 10 ...
      '9: Cloud: likely water droplets' 10 ...
      '10: Cloud: ice crystals' 10 ...
      '11: Cloud: likely ice crystals' ];

  g{'flag'}.long_definition= ...             
     ['0: No signal' 10 ...
      '1: Clean atmosphere' 10 ...
      '2: Non-typed particles/low conc.' 10 ...
      '3: Aerosol: small' 10 ...
      '4: Aerosol: large, spherical' 10 ...
      '5: Aerosol: mixture, partly non-spherical' 10 ...
      '6: Aerosol: large, non-spherical' 10 ...
      '7: Cloud: non-typed' 10 ...
      '8: Cloud: water droplets' 10 ...
      '9: Cloud: likely water droplets' 10 ...
      '10: Cloud: ice crystals' 10 ...
      '11: Cloud: likely ice crystals' ];
  
%
% Prepare for saving
%

if not_show_bad_values==1
    index = find((mask==0) | (mask==1) | (mask==2)); %find bins of no signal or clear sky
    data.ang355532(index) = NaN;        % do not save particle properties for those pixels
    data.ang5321064(index) = NaN;
    data.ang3551064(index) = NaN;
    data.pardepol(index) = NaN;
    %data.par_beta_0355(index)=NaN;
    %data.par_beta_0532(index)=NaN;
    %data.par_beta_1064(index)=NaN;
end
   

    g=store_NC_mol_beta(g,'mol_beta_0355'  , '355 nm',  molecular_beta_355,  missing_value);
    g=store_NC_mol_beta(g,'mol_beta_0532'  , '532 nm',  molecular_beta_532,  missing_value);
    g=store_NC_mol_beta(g,'mol_beta_1064'  , '1064 nm', molecular_beta_1064,  missing_value);
    g=store_NC_beta(g,'beta_0355'  , 'Quasi particle bsc at 355 nm',  data.par_beta_0355,  missing_value);
    g=store_NC_beta(g,'beta_0532'  , 'Quasi particle bsc at 532 nm',  data.par_beta_0532,  missing_value);
    g=store_NC_beta(g,'beta_1064'  , 'Quasi particle bsc at 1064 nm', data.par_beta_1064,  missing_value);
    g=store_NC_ang(g,'ang355532'   , 'Angstroem Exponent (355 - 532)' ,  data.ang355532,  missing_value);
    g=store_NC_ang(g,'ang5321064'  , 'Angstroem Exponent (532 - 1064)',  data.ang5321064,  missing_value);
    g=store_NC_ang(g,'ang3551064'  , 'Angstroem Exponent (355 - 1064)',  data.ang3551064,  missing_value);
    g=store_NC_depol(g,'particle_depolarization', 'Quasi particle depolarization ratio (532 nm)', data.pardepol, missing_value);
    g=store_NC_depol(g,'volume_depolarization', 'Calibrated volume depolarization (532 nm)', data.volume_depolarization, missing_value);
   % g=store_NC_mass_conc(g,'dust_mass_conc', 'Dust mass concetration after Poliphon', data.dust_mass_conc, missing_value);
    
 
    %%--- write global attributes to netcdf file ---%%
    g.Conventions = 'CF-1.0';
    g.system      = attribute.global.system;
    g.source      = g.system;
    g.location    = attribute.global.location;
    g.institution = attribute.global.institution;
    g.title       = attribute.global.title;
    history       = '';
    g.day   = ncshort(str2double(thedate(7:8)));
    g.month =  ncshort(str2double(thedate(5:6)));
    g.year  =  ncshort(str2double(thedate(1:4)));

    g=close(g);
   else
    disp(['No data to write for ' file_aclass]);
   end
end



    
function  g = store_NC_background(g,variable,title,data,missing_value)
    g{variable}               = ncfloat('time');
    g{variable}(:)            = data;
    g{variable}.units         = ' ' ;
    g{variable}.long_name     = [variable '(' title ')'];
    
    g{variable}.missing_value = missing_value ;
    g{variable}.FillValue_    = missing_value ;    


function  g = store_NC_beta(g,variable,title, data,missing_value)
    g{variable}               = ncfloat('time','range');
    g{variable}(:)            = data;
    g{variable}.units         = 'sr-1 Mm-1' ;
    g{variable}.units_html    = 'sr<sup>-1</sup> Mm<sup>-1</sup>' ;
    %g{variable}.long_name     = ['Quasi particle backscatter coefficient (' title ')'];
    g{variable}.long_name     = [variable ' (' title ')'];
    g{variable}.missing_value = missing_value ;
    g{variable}.FillValue_    = missing_value ;
    if strcmp(variable,'beta_1064')
      g{variable}.plot_range = [0 5e-6];%was two
    elseif  strcmp(variable,'beta_0532')
      g{variable}.plot_range = [0 10e-6];
    else
       g{variable}.plot_range = [0 1e-5]; 
    end
    g{variable}.plot_scale = 'linear';
   
    function  g = store_NC_mass_conc(g,variable,title, data,missing_value)
    g{variable}               = ncfloat('time','range');
    g{variable}(:)            = data;
    g{variable}.units         = 'mug cm⁻3' ;
    g{variable}.units_html    = 'mug<sup></sup> cm<sup>-3</sup>' ;
    %g{variable}.long_name     = ['Quasi particle backscatter coefficient (' title ')'];
    g{variable}.long_name     = [variable ' (' title ')'];
    g{variable}.missing_value = missing_value ;
    g{variable}.FillValue_    = missing_value ;
    g{variable}.plot_range = [0 400];
    g{variable}.plot_scale = 'linear';
    
    
function  g = store_NC_mol_beta(g,variable,title, data,missing_value)
    g{variable}               = ncfloat('time','range');
    g{variable}(:)            = data;
    g{variable}.units         = 'sr-1 m-1' ;
    g{variable}.units_html    = 'sr<sup>-1</sup> m<sup>-1</sup>' ;
    g{variable}.long_name     = ['Molecular backscatter coefficient (' title ')'];
    g{variable}.long_name     = [variable ' (' title ')'];
    g{variable}.missing_value = missing_value ;
    g{variable}.FillValue_    = missing_value ;
    g{variable}.plot_range = [1e-7 1e-4];
    g{variable}.plot_scale = 'logarithmic';

function  g = store_NC_depol(g,variable, title, data,missing_value)
    g{variable}               = ncfloat('time','range');
    g{variable}(:)            = data;
    g{variable}.units         = '1' ;
    g{variable}.units_html    = '1' ;
    g{variable}.long_name     = [variable ' (' title ')'];
    g{variable}.missing_value = missing_value ;
    g{variable}.FillValue_    = missing_value ;
    g{variable}.plot_range = [0 0.5];
    g{variable}.plot_scale = 'linear';
      
        
 function  g = store_NC_ang(g,variable,title, data,missing_value)
    g{variable}               = ncfloat('time','range');
    g{variable}(:)            = data;
    g{variable}.units         = '1';
    g{variable}.units_html    = '1' ;
        g{variable}.long_name     = [variable ' (' title ')'];
    g{variable}.missing_value = missing_value ;
    g{variable}.FillValue_    = missing_value ;
    g{variable}.plot_range = [0 2];
    g{variable}.plot_scale = 'linear';
     
     


function index = despeckle(data)
   [a,b] = size(data);
   index = find(~isfinite(data(1:a-2,2:b-1)) & ~isfinite(data(3:a,2:b-1)) & ~isfinite(data(2:a-1,1:b-2)) & ~isfinite(data(2:a-1,3:b)));



function param_out = remove_inner_pixels(param, index, type, value)
  if nargin<3
    type='line';
    value = NaN;
  end

  if nargin<4
    value = NaN;
  end

  switch type
    case 'line'
     param_inner = param(2:end-1,:);
     param_inner(index) = value;
     param_out = param;
     param_out(2:end-1,:) = param_inner;
    case 'despeckle'  
     param_inner = param(2:end-1,2:end-1);
     param_inner(index) = value;
     param_out = param;
     param_out(2:end-1,2:end-1) = param_inner;
  end

function droplet_bit = get_droplet_bit(height, beta)

beta0 = beta;
beta0(find(~isfinite(beta0))) = 0;

droplet_bit = zeros(size(beta));
size(beta)
dheight = median(diff(height));
nrays = size(beta,1);
ngates = length(height);

threshold_beta = 2e-5;
min_jump = 10; % Factor of 20 
jump_distance = 250; % m
jump_gates = ceil(jump_distance / dheight);
final_gate = ngates-jump_gates;
search_gates_below = ceil(100/dheight);
search_gates_above = ceil(300/dheight);
%radar_search_gates_above = ceil(2000/dheight);
diff_factor = 0.25;
%step_gates = 2;
%min_temperature = 273.16 - 40; % No droplets colder than this

%return;

disp('Calculating location of liquid water cloud droplets')
for ii = 1:nrays
  %base_liquid = Nan
  %top_liquid = Nana
   start_gate = 2;
  profile = beta0(ii,:);
  while start_gate <= final_gate
    ihighbeta = min(find(profile(start_gate:(length(profile)-search_gates_above)) > threshold_beta)) ...
	+ start_gate - 1;
    % Peak beta might be a little higher: find the highest beta in
    % this vicinity - but would this be the best approach?
    %    [] = max(profile(ihighbeta+[0:search_gates_above])
    %ihighbeta = ihighbeta + 
    if isempty(ihighbeta)
      %if (start_gate > 2)
	% Lidar extinguished - find some radar-only liquid clouds
	%top_gate = min(find(cold_bit(ii,:)));
	
    %  else
    %    if (profile(1) > threshold_beta) 
    %      droplet_bit(ii,1) = 1;
    %      top_gate = min(find(cold_bit(ii,:)));
    %      % do we need start_gate = 1 here?
    %      if (top_gate > start_gate)
    %        top_gate = max(find(~isfinite(Z(ii,start_gate:top_gate))))+start_gate-1;
    %        if ~isempty(top_gate)
    %          droplet_bit(ii, find(isfinite(Z(ii,start_gate:top_gate)))+start_gate-1) = 1;
    %        end
    %      end
    %    end
      %end
      %ii
      break
    end
    if min(profile(ihighbeta+[1:jump_gates]) ./ profile(ihighbeta)) < 1./min_jump
      % Candidate liquid water
      % Find the base
      search_start = max(1, ihighbeta-search_gates_below);
      diff_profile = diff(profile(search_start:ihighbeta));
      max_diff = max(diff_profile);
            base_liquid = min(find(diff_profile > max_diff*diff_factor))+search_start;
      % Find the top
      % First see if profile goes to zero
      top_liquid = min(find(~profile(ihighbeta+[1:search_gates_above])))+ihighbeta-1;
      if isempty(top_liquid)
        diff_profile = diff(profile(ihighbeta+[0:search_gates_above]));
        max_diff = max(-diff_profile);
        top_liquid = max(find(-diff_profile > max_diff*diff_factor))+ihighbeta-1;
      end
     droplet_bit(ii, base_liquid:top_liquid) = 1; 
    end
    start_gate = ihighbeta+1;  
      % Set to wet
    
  end
    % Previously we did this which seems very inefficient: repeated
    % analysis of the same pixels:
    %    start_gate+step_gates;
    % This seems more sensible:
    
end


%droplet_bit(find(temperature < min_temperature)) = 0; 



