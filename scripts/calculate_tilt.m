% Script to calculate max tilt of the plume by determining the max
% temperature at each depth, and finding the angle between points.
% interval can be set.
% This script uses the extracted paraview point data.
clear all


for c1=4
    for c2=1

radius1 = 5711000;  %Only look at the uppermost mantle. 5961 for umm
radius2 = 6261000;  %Cut off top to avoid looking at where plume is horizontal.
intradius = radius1 + 0.5*(radius2 - radius1);
tilt = zeros(length(1600),7);
interval= 1000;
sample = 100000;


PathStr='/home/djneuh/Desktop/dislocation_runs/m4-1_dislocation-highres/data';
PathStr2= PathStr;

%PathStr='/media/djneuh/Elements/pulsing_plumes/Thesis models/30-km/movingtop/1e6/1e6-1e6-30km-8cmtop/datap';
%PathStr2 = '/media/djneuh/Elements/pulsing_plumes/Thesis models/30-km/movingtop/1e6/1e6-1e6-30km-8cmtop/data'; %PathStr;



%PathStr=['/media/djneuh/Elements/pulsing_plumes/Thesis models/30-km/movingtop/',num2str(c1), ...
%      'e6/',num2str(c1),'e6-',num2str(c2),'e6-30km-8cmtop/datap/'];
%PathStr2=['/media/djneuh/Elements/pulsing_plumes/Thesis models/30-km/movingtop/',num2str(c1), ...
%      'e6/',num2str(c1),'e6-',num2str(c2),'e6-30km-8cmtop/data/'];
textFiles=dir([PathStr '*.csv']);


for k = 300:1600  %1000 1002 1004
  % Open point data and gather it.
  textFilename = ['point.' num2str(k) '.csv'];
  fullname = fullfile(PathStr, textFilename);
  fid=fopen(fullname, 'rt');
  Data = csvread(fullname,1);
  fclose(fid);
  
  %Open field data so we know the times.
  textFilename = ['field.' num2str(k) '.csv'];
  fullname = fullfile(PathStr2, textFilename);
  fid=fopen(fullname, 'rt');
  fdata = csvread(fullname,1);
  fclose(fid);
  %So I know how far along the script is.
  k
  
%variable to hold the max temp and position at each depth.
fp = zeros((radius2-radius1)/interval+2,4);
placement = 1;

datar= length(Data(1,:)); 
radial = zeros(length(Data),datar);
ftheta = zeros(length(Data),1);

% Convert to radial coordinates.
  for i=1:length(Data)
    radial(i,1) = atand(Data(i,datar-1)/Data(i,datar-2)); % theta
    radial(i,2) = sqrt(Data(i,datar-1)^2 + Data(i,datar-2)^2)/1000; % radius
    radial(i,3) = Data(i,15);   %nonadibatic temperature 14 or 15
  end
  
%Now get rid of any zero's.
%Only keep unique x and y points.
[C, ia, ic] = unique(radial(:,1:2), 'rows');
radial = radial(ia,:);
radial = radial(any(radial,2),:);

% Find all points along a radial depth contour and then determine
% the maximum temperature and position.
for depth = radius1:interval:radius2
placement = placement + 1;
degree_to_distance = (depth*2*pi)/360;

% If the point is within 1 km of depth, consider it part of contour.
depth = depth/1000;
radial_point = [0,0];
for i=1:length(radial)
if(round(radial(i,2)) == depth)
    radial_point(i,1) = radial(i,1); % theta
    radial_point(i,2) = radial(i,3); % temperature 
    radial_point(i,3) = depth;  % radius
end
end
radial_point = radial_point(any(radial_point,2),:);    


%If it's too early then there may be no points, if thats the case don't 
% run into an error.
%empt = isempty(Dnew);
[M, pos] = max(radial_point(:,2));

% Keep the maximum temperature along depth contour.
% Cut off at 150, anything cooler than this may not be part of the main
% plume stem and is considered an artifact.
if M > 150 
    fp(placement,1) = radial_point(pos,2);  % Max temp at depth
    fp(placement,2) = radial_point(pos,1); % distance from center
    fp(placement,3) = depth;  % radius
else
    fp(placement,1) = 0;
    fp(placement,2) = 0;
    fp(placement,3) = 0;
end

clearvars -except radial datar PathStr radius1 radius2 p plume g a Data k theta fp frame F tilt fdata PathStr2 v interval placement sample c1 c2

end

% Now we have a temp and position for each depth.
fp = fp(any(fp,2),:);
    if(k==0)
      time=0;
    else 
      %convert to mya)
      time=fdata(2)/(1e6);
    end

% If we have values (i.e. the plume is in the upper mantle), and there
% enough values to indicate the plume then do this.
empt2 = isempty(fp);
if length(fp) < 15
    empt2 = 1;
end

if empt2 == 0
% check for outliers.
for i = 1:(length(fp(:,1))-1)
    % If we have sudden 60 km increase over 5 km, it's likely noise.
    % For now set only theta to 0 to not mess up difference in subsequent
    % iterations.
    val = fp(i,4);
    
    difference(i) = abs(fp(i+1,2)-fp(i-val,2))*(fp(i-val,3)*2*pi)/360;
    if( difference(i) > 60)
          fp(i+1,4) = fp(i,4)+1;
    else
          fp(i+1,4) = 0;
    end
end
%stop
% remove outliers.
for i = 1:(length(fp(:,1)))   
    %Set depths to zero now.
    if( fp(i,4) > 0)
          fp(i,1) = 0;
          fp(i,2) = 0;
          fp(i,3) = 0;
          fp(i,4) = 0;
     end
end
fp = fp(any(fp,2),:);

if(length(fp(:,1)) < 8)
    tilt(k,4) = 0; % good one
    continue
end

% Find the distance from 45 degrees, based on the radius.
fp(:,5) = abs(45 - fp(:,2));
fp(:,6) = fp(:,5).*fp(:,3)*2*pi/360;

% Fit a function to our points.
inter = 5;  % km
 x = min(fp(:,3)):inter:max(fp(:,3));   
[p,~,mu] = polyfit(fp(:,3),fp(:,6), 8);
y_fit = polyval(p,x,[],mu);


for i=1:length(y_fit)-1  
    % Find distance between points, because it's from 45 degrees
    % distance should always be positive.
    dd = abs(y_fit(i+1) - y_fit(i));
    
    % Convert to 
    angle(i) = 90 - atand(inter/dd);
end

for i=1:length(fp(:,1))-1  
    % Find distance between points, because it's from 45 degrees
    % distance should always be positive.
    dd2 = abs(fp(i+1) - fp(i));
    
    % Convert to 
    angle2(i) = 90 - atand(inter/dd2);
end

% Remove top and bottom points as they often are not robust.
[a1, pos1] = max(angle); 
[a2, pos2] = max(angle2); 




avg_interval = 25;
avg_length = round(avg_interval/inter);
pl = 1;
for i=1:length(angle)
    if(i>=pos1-2 && i<=pos1+2)
        avg(pl) = angle(i);
        pl = pl + 1;
    end        
end

for i=1:length(angle2)
    if(i>=pos1-2 && i<=pos1+2)
        avg3(pl) = angle2(i);
        pl = pl + 1;
    end        
end

for i=3:length(angle)-2
        avg2(i) = (angle(i-2)+ angle(i-1)+ angle(i)+ angle(i+1)+ angle(i+2))/5;
end

for i=3:length(angle2)-2
        avg4(i) = (angle2(i-2)+ angle2(i-1)+ angle2(i)+ angle2(i+1)+ angle2(i+2))/5;
end

tilt(k,1) = k;
tilt(k,2) = time;
tilt(k,3) = mean(avg);
tilt(k,4) = max(avg2);

[a2, pos2] = max(avg2); 
tilt(k,5) = 5711+pos2*5;

%tilt(k,5) = max(avg4);



% %%%%%%%%%%%Plot our function and how it fits points.
% clf(figure(1))
% figure(1)
% axis equal
% plot(y_fit,x,'k','LineWidth',1.5)
% hold on
% scatter(y_fit(pos1+2),x(pos1+2),'filled')
% scatter(fp(:,2),fp(:,3))
% xlabel('Distance from center (km)')
% ylabel('Radius (km)')
% pause(0.1)
%max_angle = a1

%%%%%%%%%%%%Plot showing all point and then the maximum temp line we found.
 %theta = 45-atand(x(pos1+2)/y_fit(pos1+2));
 % Find the distance from 45 degrees, based on the radius.
 %radial(:,5) = abs(45.0 - radial(:,1));
 %radial(:,6) = radial(:,5).*radial(:,2)*2*pi/360;
 
%  for i =1:length(radial)
%  if(radial(i,3) > 150)
%      pf(i,1) = radial(i,1);
%      pf(i,2) = radial(i,2);
%      pf(i,3) = radial(i,3);
%  end
%  end
%  [a1, pos1] = max(avg2); 
%  pf(:,5) = abs(45.0 - pf(:,1));
%  pf(:,6) = pf(:,5).*pf(:,2)*2*pi/360;
%  
%  x2 = [0 1400];
%  y2 = [x(pos1+2) x(pos1+2)];
%  yt = [6336 6336];
%  y410 = [5961 5961];
%  clf(figure(2))
%  figure(2)
%  scatter(pf(:,6),pf(:,2),100,pf(:,3),'.')
%  hold on
%  %scatter(fp(:,6),fp(:,3),'k.')
%  plot(y_fit,x,'k','LineWidth',1.5)
%  %set(gca, 'XDir','reverse')
%  plot(x2,y2,'r--','LineWidth',1.5)
%  plot(x2,yt,'k','LineWidth',1.5)
%  plot(x2,y410,'k:','LineWidth',1)
%  axis equal
%  %text(50,6100,['Max tilt: ', num2str(round(max(avg2)))],'Color','black','FontSize',16,'FontWeight','bold')
%  xlabel('Distance (km)')
%  ylabel('Radius (km)')
%  text(50,6375,['Top boundary'],'FontSize',10)
%  text(1100,6000,['410 transition'],'FontSize',10)
%  xlim([0 1400]);
%  ylim([5900 6400]);
%  %stop
%  
%  pause(0.1)
 end

%clearvars -except radius1 radius2 k g a PathStr mean_Ra F frame tilt  PathStr2 v interval
end

tilt = tilt(any(tilt,2),:);
save(['tilt_2021-hr-',num2str(c1),'-',num2str(c2)], 'tilt');
%clearvars -except c1 c2 tilt
clear
    end
end



