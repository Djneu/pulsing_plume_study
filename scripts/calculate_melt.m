%%Script for calculating melt in m^3/m*s%%
% For this script load the statistics file and volume flux .mat
clearvars -except statistics final410


%Set max year as the year max volume flux occurs
[M,I] = max(final410);
maxyear = final410(I(2),1)
interval = 25;
low = 70 - interval/2; %100  % I should check what number I used for periodicities.
high = 200 + interval/2;  % High number so we see everything that happens
year_to_second = 365*24*60*60;

%Make new matrix that only keeps data between a given low and high around the
%max volume flux. (Generally max is when plume first
%passes the chosen depth.)
for i = 1:length(final410)
    if final410(i,1) >= (maxyear+ low) && final410(i,1) <= (maxyear+ high)
        time200(i,1) = final410(i,1);
        time200(i,2) = final410(i,2);
    else
        time200(i,:) = 0;
    end
end
%Get rid of all the zero elements.
time200 = time200(any(time200,2),:);

%Set first element as 0 so we can track over the total 205 Myr of data we
%kept.
time200(:,1) = time200(:,1)-time200(1,1);


%Find dF and set correct years to look at
for i=2:length(statistics)
    
    %Subtract melt from previous timestep so we know how much melt is
    %produced in the given timestep, and then convert to m^3/ms. statistics 3
    %is the timestep size.
    meltt(i) = (statistics(i,16) - statistics(i-1,16))/(statistics(i,3)*year_to_second);
    
    % Now, if the statistics is within our given time range, add it to
    % variable we'll use to calculate melt from a pulse. Times for this if
    % are given in myr.
    if statistics(i,2)/1e6 >= (maxyear+low) && statistics(i,2)/1e6 <= (maxyear+high)      
        meltf(i,1) = statistics(i,2);  %timestep
        meltf(i,2) = meltt(i);         %melt produced in m^3/ms
        
        %If we want to know what the original melt looked like.
        melt_orig(i,1) = statistics(i,16);
    end
end

%get rid of zeroes and shift time to zero.    
meltf = meltf(any(meltf,2),:);
meltf(:,1) = meltf(:,1)-meltf(1,1);

scatter(meltf(:,1)/1e6,meltf(:,2),'.')
%hold on

meltf(:,2) = medfilt1(meltf(:,2),10);   %%generally 7
%scatter(meltf(:,1)/1e6,meltf(:,2),'.')
hold on

%Shift up everything so we don't deal with negatives.
%meltf(:,2) = meltf(:,2);

%Find the minimum's between pulses, we'll then use this to determine where
%a pulse starts and ends, giving a time range to integrate for each pulse.
%TF = islocalmin(meltf(:,2),'MinSeparation',3e2,'MinProminence', 2e-6); %'MinProminence', 5*10^8, 'MinSeparation',minutes(45),'SamplePoints',t  'MaxNumExtrema', 5
TF = islocalmin(meltf(:,2),'MinSeparation',2e2,'MinProminence', 8e-7);
%[val, TF] = findpeaks(-meltf(:,2),'MinPeakDistance', 400,'Threshold',1e-9);
meltmin(:,1) = meltf(TF,1);
meltmin(:,2) = meltf(TF,2);

scatter(meltmin(:,1)/1e6,meltmin(:,2),'ro','filled')
hold on
scatter(meltf(:,1)/1e6,meltf(:,2),'k.')




% Now we'll integrate each pulse by going from minimum to minimum. This can
% give us the total melt produced for each pulse.
for i=2:length(meltmin(:,1))
% Get an equation for a line between the two minimums to isolate a pulse.

[coefficients c] = polyfit([meltmin(i-1,1), meltmin(i,1)], [meltmin(i-1,2), meltmin(i,2)], 1);
a = coefficients (1);
b = coefficients (2);

% a = (meltmin(i,2) - meltmin(i-1,2))/(meltmin(i,1) - meltmin(i-1,1));
% 
% if(meltmin(i-1,1) > meltmin(i,1))
%    b = meltmin(i,2);
% else
%    b = meltmin(i-1,2);
% end

%stop




    
    for j=1:length(meltf(:,1))
               
        if(meltf(j,1) > meltmin(i-1,1) && meltf(j,1) < meltmin(i,1))
            F(j) = meltf(j,2)- (meltf(j,1)*a+b);
            dt(j) = (meltf(j+1,1) - meltf(j-1,1))/2;
            t(j) = meltf(j,1);
            yy(j) = (meltf(j,1)*a+b);
        elseif meltf(j,1) == meltmin(i,1)
            F(j) = meltf(j,2)- (meltf(j,1)*a+b);
            dt(j) = (meltf(j,1)- meltf(j-1,1))/2;
            t(j) = meltf(j,1);
            yy(j) = (meltf(j,1)*a+b);
        elseif meltf(j,1) == meltmin(i-1,1)
            F(j) = meltf(j,2)- (meltf(j,1)*a+b);
            dt(j) = (meltf(j+1,1) - meltf(j,1))/2;
            t(j) = meltf(j,1);
            yy(j) = (meltf(j,1)*a+b);
        end      
    end
    
    tt = t;
    tt = tt(:, any(tt,1));
    yy = yy(:, any(yy,1));
    
%      for j=1:length(F)
%          if(F(j) < 0)
%             F(j) = 0;
%         end
%      end
    
    
    %Plot area of each pulse
    %min(dt);
    %figure(1)
    xlabel('Time (myr)');
    ylabel('F (m^3/ms)');
    title('Melt Flux (410-km = 4 MPa/K @ 4 cm/yr) ');
%    scatter(t/1e6,F,'.');
    plot(tt/1e6,yy,'r');
   hold on
% 
%     if mod(i,2)==0
%      area(t/1e6,F,'FaceAlpha',1,'EdgeAlpha', 0,'FaceColor',[0.6 0.8 1])
%      scatter(t/1e6,F,'k.');
%     else
%      area(t/1e6,F,'FaceAlpha',1,'EdgeAlpha', 0,'FaceColor',[0.2 0.6 0.8])
%      scatter(t/1e6,F,'.','MarkerEdgeColor',[0.5 0.5 0.5]);
%     end
%     
    
    
    %if i==3
    %    stop
    %end
   
    
    % Multiply melt produced per time by dt to find total melt of the
    % pulse.
    for jj=1:length(dt)
        local_integration(jj) = F(jj)*dt(jj)*year_to_second;
    end
    
    %remove zeroes for averaging.
    %F = F(:, any(F,1));
    %dt = dt(:, any(dt,1));
    
    
   % Our averaged pulse size in km^3/km
   melt_integration(i-1) = sum(local_integration)*1e-6; %*1e-6; if correcting to km^3/km
   flux_integration(i-1) = sum(local_integration)/((meltmin(i,1)-meltmin(i-1,1))*year_to_second);
   
   % And our average melt production in m^3/ms
   melt_average(i-1) = mean(F);
   melt_peak(i-1) = max(F);

   clear local_integration F t dt index yy
    
end

pav = mean(melt_integration);  %%in km^3/km
for ji = 1:length(melt_integration)
  melt_dif(ji) = pav - melt_integration(ji);
end
%stop
% Now put it all into a matrix I can copy.
 final(1) = round(pav,2)
 final(3) = round(min(melt_integration),2);
 final(4) = round(max(melt_integration),2);
 final(2)  = round(sqrt( mean( (melt_dif).^2 )),2 );

 final(5) = mean(flux_integration)  %mean(melt_average);
 %final(6) = min(flux_integration);
 %final(7) = max(flux_integration)
 
 %final(8) = mean(melt_peak);
 %final(9) = min(melt_peak);
 %final(10) = max(melt_peak)

clearvars -except final



