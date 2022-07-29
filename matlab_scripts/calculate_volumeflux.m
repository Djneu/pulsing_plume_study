%Script to load data from aspect, and then find the volume flux for each
%timestep.

clear
tic 

% Set a depth contour, and then use 40 km of data surrounding that contour
% for interpolating.
intdepth = 175; % in km
intradius = (6336-intdepth)*1e3;
radius1 = intradius - 20e3; %6096000;  %lower end of radius 5966000   6096000   6291000   !6216
radius2 = intradius + 20e3; %6136000;  %high end of radius 6026000    6136000   6331000   !6256

%Find the distance a degree is at given radius
degree_to_distance = (intradius*2*pi)/360;
final410 = zeros(1600,2);

% Set range (in theta) to consider volume flux around maximum nonadiabatic temperature
thetarange = 1;

PathStr=['/home/djneuh/final_runs/r2tz31/data/'];
name = "vr_r2_tz31_v8_t550_p1"

%create a variable linked to the directory with data files
textFiles=dir([PathStr '*.csv']);

%Here we run through all output paraview files to get the time (field) data
%and all point data. For this we always save field data as field, and point
%data as point
for k = 1:(floor(length(textFiles)/2)-1) 
  final_temp = zeros(1,2);  

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
  
  %Do the same for field data.
  textFilename = ['field_' num2str(k) '.csv'];
  fullname = fullfile(PathStr, textFilename);
  fid=fopen(fullname, 'rt');
  fdata = csvread(fullname,1);
  fclose(fid);
  
  %Outputting k because I'm too impatient to wait for the program to fully
  %run without knowing where it is.
  k
  
 
  
  %Setting up a matrix with 0's for run speed, and then keeping only the
  %points that are in our chosen radius range.
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
    Dd(i,1) = test2(i,1);         %%1
    Dd(i,2) = test2(i,2);         %%2
    Dd(i,3) = test2(i,datar-2);        %%3
    Dd(i,4) = test2(i,datar-1);        %%4
end

%Get rid of any duplicate points
Dd = unique(Dd, 'rows');

%Calculate velocity in radial direction for each point
theta = zeros(length(Dd),1);
Vr = zeros(length(Dd),1);

for j=1:length(Dd)
  theta(j,:)=atand(Dd(j,4)/Dd(j,3));
  Vr(j,:)=Dd(j,1)*cosd(theta(j))+Dd(j,2)*sind(theta(j));
end
 
%Get rid of any duplicates again.
theta = unique(theta);

%calculate our x and y for the radius we want.
xq = zeros(length(theta),1);
yq = zeros(length(theta),1);

for n=1:length(theta)
    xq(n,:) = intradius*cosd(theta(n)); 
    yq(n,:) = intradius*sind(theta(n));
end

%interpolate radial velocities at chosen radius
Vrin = zeros(length(xq),3);
Vrin(:,1) = griddata(Dd(:,3), Dd(:,4), Vr, xq, yq);
Vrin(:,2) = theta;
Vrin(any(isnan(Vrin), 2), :) = []; %Get rid of any NaN points

% Find the maximum nonadiabatic temperature along chosen 
% contour to consider as plume centr.
[M, pos] = max(Vrin(:,1));
ptheta = Vrin(pos,2);

% Apply theta range.
mantle_visc = zeros(length(Vrin),2);
for df=1:length(Vrin(:,1))
    if Vrin(df,2) < ptheta+thetarange && Vrin(df,2) > ptheta-thetarange
        mantle_visc(df,1) = Vrin(df,1);
        mantle_visc(df,2) = Vrin(df,2);
    end
end

%We're only interested in outward velocity, so we only keep positive radial
%velocities.
Vfinal = zeros(length(Vrin),3);
for j=1:length(Vrin)
    % We don't want to include negative (inward) velocities.
    if(Vrin(j,1)<0)
      Vfinal(j,1) = 0; %velocity
      Vfinal(j,2) = mantle_visc(j,2);  %Theta value
    else
      Vfinal(j,1) = mantle_visc(j,1);
      Vfinal(j,2) = mantle_visc(j,2);
    end
end

%If we have positive velocites in this timestep set to zero
% so it isn't included.
ex = exist('Vfinal', 'var');
if ex==1
    Vfinal(Vfinal(:,2)==0, :)= [];
else 
    %if we don't have positive, make a zero matrix so we don't get an
    %error when we try to use it.
    Vfinal=[0,0;0,0];
end

%Now find the change in theta between each point.
dtheta = zeros(length(Vfinal),1);
for p=1:length(Vfinal)
    if p==1
        dtheta(p,:) = (Vfinal(p+1,2)-Vfinal(p,2))/2;
    elseif p == length(Vfinal)
        dtheta(p, :) = (Vfinal(p,2) - Vfinal(p-1,2))/2;
    else
        dtheta(p,:) = (Vfinal(p+1,2) - Vfinal(p-1,2))/2;
    end
end


%Now calculate the volume flux for each point.
ff = zeros(length(Vfinal),1);
for n=1:length(Vfinal)
    ff(n,:) = Vfinal(n,1)*dtheta(n)*degree_to_distance;
end

%If we have a volume flux, put it into final matrix as a sum of all the
%points.
if(sum(ff)>0)
    final_temp(1,2) = sum(ff);
else
    final_temp(1,2) = 0;
end


%set first timestep to 0 since output data is weird.
  if(k==0)
    final_temp(1,1)=0;
  else 
    %convert to mya)
    final_temp(1,1)=fdata(2)/(1e6);
  end

  final410(k+1,:) = final_temp;

%Clear all the variables that may give us trouble, leave any that we need.
%clearvars -except final410 PathStr intradius degree_to_distance radius1 radius2 final660 btheta fffs x410 x660 v PathStrf PathStrp
 
end

final410 = final410(any(final410,2),:);
save(name, 'final410')

%Short plot just to test and make sure things are working.
%figure(1),
%plot(final410(:,1),final410(:,2), 'LineWidth', 1.5)
%title('Time vs. Volume flux 3e6-2e6 (depth 350 km) ');
%xlabel('Time (Myr)');
%ylabel('Volume flux (m^3/yr)');

toc