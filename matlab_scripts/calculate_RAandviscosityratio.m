% This scripts uses the extracted paraview point data.

clear
tic

% Set a depth contour, and then use 40 km of data surrounding that contour
% for interpolating.
intdepth = 340; % in km
intradius = (6336-intdepth)*1e3;
radius1 = intradius - 20e3; %6096000;  %lower end of radius 5966000   6096000   6291000   !6216
radius2 = intradius + 20e3; %6136000;

%Find the distance a degree is at given radius
degree_to_distance = (intradius*2*pi)/360;
final410 = zeros(1600,12);


PathStr=['/home/djneuh/final_runs/r2tz31/data/'];
name = "vr_r2_tz31_v8_t550_p1"

%create a variable linked to the directory with data files
textFiles=dir([PathStr '*.csv']);

%Here we run through all output paraview files to get the time (field) data
%and all point data. For this we always save field data as field, and point
%data as point
for k = 300:(floor(length(textFiles)/2)-1)  
  final_temp = zeros(1,12);

  %We open the file with given name and number, and then assign its path
  %and name to variable fullname.
  textFilename = ['point_' num2str(k) '.csv'];
  fullname = fullfile(PathStr, textFilename);
  
  %Here we open the text file as read only (rt), and then read the info
  %into Data.
  fid=fopen(fullname, 'rt');
  Data = csvread(fullname,1);
  
  %Close the text file
  fclose(fid);
  
  %Do the same for field (time) data.
  textFilename = ['field_' num2str(k) '.csv'];
  fullname = fullfile(PathStr, textFilename);
  fid=fopen(fullname, 'rt');
  fdata = csvread(fullname,1);
  fclose(fid);
  
  %Outputting k because I'm too impatient to wait for the program to fully
  %run without knowing where it is.
  k
  
  % keeping only the points that are in our chosen radius range.
  datar= length(Data(1,:));
  r = zeros(length(Data),1);
  test2 = zeros(length(Data),datar);           
  for i=1:length(Data)
    r(i,:)=sqrt(Data(i,datar-2)^2+Data(i,datar-1)^2);           %%23 and 24 (original) 27 and 28 (second runs) 25 and 26 (new runs)
   if(r(i,:)>radius1 && r(i,:)<radius2)                       
       test2(i,:)=Data(i,:);
    end
  end
  
  
%Now get rid of any zero's
test2 = test2(any(test2,2),:);
Dd = zeros(length(test2),4);

%Take only needed point info and put into new matrix, make sure to do this
%before so we don't mix up point info.
for i=1:length(test2)
    Dd(i,1) = test2(i,5);     %temperature
    Dd(i,2) = test2(i,7);    %viscosity, 10 if melt is included, otherwise 9
    Dd(i,3) = test2(i,datar-2);        %%3 velocity x
    Dd(i,4) = test2(i,datar-1);        %%4 velocity y
    
    Dd(i,5) = test2(i,6);    %density   normall 7
    Dd(i,6) = test2(i,8);   %ad_density  norm 13
    Dd(i,7) = test2(i,4);  %pressure
    Dd(i,8) = test2(i,9);  %nonad_temp
    Dd(i,9) = test2(i,10);  %nonad_temp
    
end

%Get rid of any duplicate points
Dd = unique(Dd, 'rows');

%Calculate velocity in radial direction for each point
theta = zeros(length(Dd),1);

for j=1:length(Dd)
  theta(j,:)=atand(Dd(j,4)/Dd(j,3));
end
 
%Get rid of any duplicates again.
theta = unique(theta);

%calculate our x and y for the depth contour we want.
xq = zeros(length(theta),1);
yq = zeros(length(theta),1);

for n=1:length(theta)
    xq(n,:) = intradius*cosd(theta(n)); 
    yq(n,:) = intradius*sind(theta(n));
end

%interpolate radial velocities at chosen radius
Vrin = zeros(length(xq),3);
Vrin(:,1) = griddata(Dd(:,3), Dd(:,4), Dd(:,1), xq, yq);   %Temperature
Vrin(:,2) = theta;
Vrin(:,3) = griddata(Dd(:,3), Dd(:,4), Dd(:,2), xq, yq);   %viscosity
Vrin(:,4) = griddata(Dd(:,3), Dd(:,4), Dd(:,5), xq, yq);   %density
Vrin(:,5) = griddata(Dd(:,3), Dd(:,4), Dd(:,6), xq, yq);   %ad density
Vrin(:,6) = griddata(Dd(:,3), Dd(:,4), Dd(:,7), xq, yq);   %pressure
Vrin(:,7) = griddata(Dd(:,3), Dd(:,4), Dd(:,8), xq, yq);   %nonad temp
Vrin(:,8) = griddata(Dd(:,3), Dd(:,4), Dd(:,9), xq, yq);   %diffusion
Vrin(any(isnan(Vrin), 2), :) = []; %Get rid of any NaN points

[p_temp,I] = max(Vrin(:,1));
[pna_temp,I2] = max(Vrin(:,7));

ptheta = Vrin(I,2);
visc_plume = Vrin(I,3);

% Set the distance interval to integrate mantle properties.
mantle_visc = zeros(length(Vrin),6);
for df=1:length(Vrin(:,1))
    if (Vrin(df,2) - ptheta) > 2 && (Vrin(df,2) - ptheta) < 12
        mantle_visc(df,1) = Vrin(df,3);
        mantle_visc(df,2) = Vrin(df,2);
        mantle_visc(df,3) = Vrin(df,4);
        mantle_visc(df,4) = Vrin(df,6); %pressure
        mantle_visc(df,5) = Vrin(df,1);  %temp
        mantle_visc(df,6) = Vrin(df,7);  %conductivity non_ad
        mantle_visc(df,7) = Vrin(df,8);  %diffusion
    end
end

mantle_visc = mantle_visc(any(mantle_visc,2),:);

% If we need conductivity calculate it since it wasn't output.
% Not used in study.
c0 = 2.47;
c1 = 0.33;
c2 = 0.48;
conduc = zeros(length(mantle_visc(:,1)),1);
for cc = 1:length(mantle_visc(:,1))
  conduc(cc,1) = (c0+c1*mantle_visc(cc,4)*1e-9)*((300/mantle_visc(cc,5))^c2);
end

% Find change in theta between points
dtheta = zeros(length(mantle_visc(:,1)),1);
for p=1:length(mantle_visc(:,1))
    if length(mantle_visc(:,1)) == 1
        dtheta(p,:) = 0;
    elseif p==1 
        dtheta(p,:) = (mantle_visc(p+1,2)-mantle_visc(p,2))/2;
    elseif p == length(mantle_visc(:,1))
        dtheta(p,:) = (mantle_visc(p,2) - mantle_visc(p-1,2))/2;
    else
        dtheta(p,:) = (mantle_visc(p+1,2) - mantle_visc(p-1,2))/2;
    end
end

% Integrate along interval.
ff = zeros(length(mantle_visc(:,1)),1);
ff2 = zeros(length(mantle_visc(:,1)),1);  
ff3 = zeros(length(mantle_visc(:,1)),1);   
ff4 = zeros(length(mantle_visc(:,1)),1);  
ff5 = zeros(length(mantle_visc(:,1)),1);   
ff6 = zeros(length(mantle_visc(:,1)),1);   
for n=1:length(mantle_visc(:,1))
    ff(n,:) = mantle_visc(n,1)*dtheta(n)*degree_to_distance;
    ff2(n,:) = mantle_visc(n,3)*dtheta(n)*degree_to_distance;
    ff3(n,:) = mantle_visc(n,6)*dtheta(n)*degree_to_distance;
    ff4(n,:) = conduc(n,1)*dtheta(n)*degree_to_distance;
    ff5(n,:) = mantle_visc(n,5)*dtheta(n)*degree_to_distance;
    ff6(n,:) = mantle_visc(n,7)*dtheta(n)*degree_to_distance;
end

% Find average for each property.
if(sum(ff)>0)
    visc_mantle = sum(ff)/(sum(dtheta)*degree_to_distance);
else
    visc_mantle = 0;
end

if(sum(ff2)>0)
    density = sum(ff2)/(sum(dtheta)*degree_to_distance);
else
    density = 0;
end

if(sum(ff3)>0)
    conductivity = sum(ff3)/(sum(dtheta)*degree_to_distance);
else
    conductivity = 0;
end

if(sum(ff4)>0)
    conductivity_calculated = sum(ff4)/(sum(dtheta)*degree_to_distance);
else
    conductivity_calculated = 0;
end

if(sum(ff5)>0)
    mantle_temp = sum(ff5)/(sum(dtheta)*degree_to_distance);
else
    mantle_temp = 0;
end

if(sum(ff6)>0)
    diffusion = sum(ff6)/(sum(dtheta)*degree_to_distance);
else
    diffusion = 0;
end


  visc_ratio = visc_mantle/visc_plume;
  meanrho = abs(Vrin(I,5) - Vrin(I,4));
  kk = conductivity_calculated/(density*1250);   %thermal diffusivity from average
  g = 10;
  a = 25*1000;
  Ra = (meanrho*g*a^3)/(visc_mantle*kk);

  final_temp(1,2) = visc_ratio;
  final_temp(1,3) = visc_mantle;
  final_temp(1,4) = visc_plume;
  final_temp(1,5) = p_temp;
  final_temp(1,6) = ptheta;
  final_temp(1,7) = Ra;
  final_temp(1,10) = mantle_temp;
  final_temp(1,11) = pna_temp;
  final_temp(1,12) = diffusion;
  
  final_temp(1,9) = meanrho;
  
%set first timestep to 0 since output data is weird.
  if(k==0)
    final_temp(1,1)=0;
  else 
    %convert to mya)
    final_temp(1,1)=fdata(2)/(1e6);
  end
  
final410(k+1,:) = final_temp;
%Clear all the variables that may give us trouble, leave any that we need.
%clearvars -except final410 PathStr intradius degree_to_distance radius1 radius2 final660 btheta fffs x410 x660 PathStr2

end

final410 = final410(any(final410,2),:);
save(name, 'final410')

toc