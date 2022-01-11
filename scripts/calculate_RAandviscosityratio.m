% This scripts uses the extracted paraview point data.

%Script to load data from aspect, and then find the volume flux for each
%timestep.
for x410 = 4:4
    for x660 = 1:1

%Set radius range for gathering points, and find the radius we will
%interpret at by being midway between.      
radius1 =6016000;  %lower end of radius 5966000   6096000   6216000
radius2 = 6056000;  %high end of radius 6026000    6136000   6256000
intradius = radius1 + 0.5*(radius2 - radius1);
%Find the distance a degree is at given radius
degree_to_distance = (intradius*2*pi)/360;
final410 = zeros(1600,2);



%go into directory with data
%/media/derek/Seagate Expansion Drive/Seagate/Thesis models/30-km/movingtop/2e6/2e6-2e6-30km-8cmtop/data
%PathStr=['/home/djneuh/Desktop/dislocation_runs/m' num2str(v) '-1_dislocation-8cm/data'];
PathStr=['/media/djneuh/Elements/pulsing_plumes/Thesis models/30-km/movingtop/',num2str(x410), ...
      'e6/',num2str(x410),'e6-',num2str(x660),'e6-30km-6cmtop/datap/'];
PathStr2=['/media/djneuh/Elements/pulsing_plumes/Thesis models/30-km/movingtop/',num2str(x410), ...
      'e6/',num2str(x410),'e6-',num2str(x660),'e6-30km-6cmtop/data/'];

%create a variable linked to the directory with data files
textFiles=dir([PathStr '*.csv']);

%Here we run through all output paraview files to get the time (field) data
%and all point data. For this we always save field data as field, and point
%data as point
for k = 300:(length(textFiles)-1)
    
  %We open the file with given name and number, and then assign its path
  %and name to variable fullname.
  textFilename = ['point.' num2str(k) '.csv'];
  fullname = fullfile(PathStr, textFilename);
  
  %Here we open the text file as read only (rt), and then read the info
  %into Data.
  fid=fopen(fullname, 'rt');
  Data = csvread(fullname,1);
  
  %Close the text file
  fclose(fid);
  
  %Do the same for field data.
  textFilename = ['field.' num2str(k) '.csv'];
  fullname = fullfile(PathStr2, textFilename);
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
    Dd(i,1) = test2(i,5);     %temperature
    Dd(i,2) = test2(i,10);    %viscosity, 10 if melt is included, otherwise 9
    Dd(i,3) = test2(i,datar-2);        %%3
    Dd(i,4) = test2(i,datar-1);        %%4
    
    Dd(i,5) = test2(i,7);    %density   normall 7
    Dd(i,6) = test2(i,13);   %ad_density  norm 13
    Dd(i,7) = test2(i,4);  %pressure
    Dd(i,8) = test2(i,16);  %conductivity
    
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

%calculate our x and y for the radius we want.
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
Vrin(:,7) = griddata(Dd(:,3), Dd(:,4), Dd(:,8), xq, yq);   %conductivity (if there is any...)
Vrin(any(isnan(Vrin), 2), :) = []; %Get rid of any NaN points

[p_temp,I] = max(Vrin(:,1));
ptheta = Vrin(I,2);
visc_plume = Vrin(I,3);


%%old method to get mantle viscosity
% for df=1:length(Vrin(:,1))
%     if df == length(Vrin(:,1))
%         break
%     end
%     if Vrin(df,2) - ptheta > 6
%         mtheta = Vrin(df,2);
%         break
%     end
% end

mantle_visc = zeros(length(Vrin),6);
for df=1:length(Vrin(:,1))
    if (Vrin(df,2) - ptheta) > 2 && (Vrin(df,2) - ptheta) < 12
        mantle_visc(df,1) = Vrin(df,3);
        mantle_visc(df,2) = Vrin(df,2);
        mantle_visc(df,3) = Vrin(df,4);
        mantle_visc(df,4) = Vrin(df,6); %pressure
        mantle_visc(df,5) = Vrin(df,1);  %temp
        mantle_visc(df,6) = Vrin(df,7);  %conductivity
    end
end

mantle_visc = mantle_visc(any(mantle_visc,2),:);

c0 = 2.47;
c1 = 0.33;
c2 = 0.48;
%(c0[ol_index]+(c1[ol_index]*adiabatic_pressure*1e-9))*pow((300/temperature),c2[ol_index]);
conduc = zeros(length(mantle_visc(:,1)),1);
for cc = 1:length(mantle_visc(:,1))
  conduc(cc,1) = (c0+c1*mantle_visc(cc,4)*1e-9)*((300/mantle_visc(cc,5))^c2);
end

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

ff = zeros(length(mantle_visc(:,1)),1);
ff2 = zeros(length(mantle_visc(:,1)),1);   %density
ff3 = zeros(length(mantle_visc(:,1)),1);   %conductivity
ff4 = zeros(length(mantle_visc(:,1)),1);   %conductivity
for n=1:length(mantle_visc(:,1))
    ff(n,:) = mantle_visc(n,1)*dtheta(n)*degree_to_distance;
    ff2(n,:) = mantle_visc(n,3)*dtheta(n)*degree_to_distance;
    ff3(n,:) = mantle_visc(n,6)*dtheta(n)*degree_to_distance;
    ff4(n,:) = conduc(n,1)*dtheta(n)*degree_to_distance;
end

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
%if df ~= length(Vrin(:,1))
  %m_temp = Vrin(df,1);
  %visc_mantle = Vrin(df,3);
  visc_ratio = visc_mantle/visc_plume;
  meanrho = abs(Vrin(I,5) - Vrin(I,4));
  kk = conductivity_calculated/(density*1250);   %thermal diffusivity from average
  g = 10;
  a = 25*1000;
  Ra = (meanrho*g*a^3)/(visc_mantle*kk);
%end

  final410(k+1,2) = visc_ratio;
  final410(k+1,3) = visc_mantle;
  final410(k+1,4) = visc_plume;
  final410(k+1,5) = p_temp;
  final410(k+1,6) = ptheta;
  final410(k+1,7) = Ra;
  
  final410(k+1,8) = conductivity_calculated;
  %final410(k+1,9) = conductivity;
  final410(k+1,9) = meanrho;
  
%set first timestep to 0 since output data is weird.
  if(k==0)
    final410(k+1,1)=0;
  else 
    %convert to mya)
    final410(k+1,1)=fdata(2)/(1e6);
  end
  
%Clear all the variables that may give us trouble, leave any that we need.
clearvars -except final410 PathStr intradius degree_to_distance radius1 radius2 final660 btheta fffs x410 x660 PathStr2


  
end


final410 = final410(any(final410,2),:);

%save(['ra_300km_',num2str(x410),'-',num2str(x660)], 'final410');
%save([num2str(v) '_ra_new'], 'final410')
%clearvars -except x410 x660
    end
end

%Short plot just to test and make sure things are working.
%figure(1),
%plot(final410(:,1),final410(:,2), 'LineWidth', 1.5)
%title('Time vs. Volume flux 3e6-2e6 (depth 350 km) ');
%xlabel('Time (Myr)');
%ylabel('Volume flux (m^3/yr)');
