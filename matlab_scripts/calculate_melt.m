%%Script for calculating melt%%
% For this script load the statistics file and volume flux .mat
clearvars -except statistics final410
clf
name = '_r2_tz31_v8_t550_p1.mat';
umpt = 202;   %UMP initiation time

load(['/home/djneuh/Desktop/movies/final_runs/data_files/volume_flux/newdata/vf', name])
load(['/home/djneuh/Desktop/movies/final_runs/data_files/statistics/stats', name])

%Set max year as the year max volume flux occurs
[M,I] = max(final410);
maxyear = final410(I(2),1)
interval = 25;
year_to_second = 365*24*60*60;


% Both of these runs do not have a max flux that coincides with plume
% arrival, so we adjust here.
if contains('vf_r2_tz31_v8_t550_p10.mat',name2)
    maxyear = 140;
end

if contains('vf_r2_tz11_v8_t550_p1.mat',name2)
    maxyear = 118.3;
end

low =  umpt - maxyear; 
high = 400-maxyear-interval/2; 


% Smooth the statistics integrated global melt data.
statistics_smoothed = smoothdata(statistics(:,16),'movmean',10);

%Find dF and set correct years to look at
for i=2:length(statistics-1)   
    % Now, if the statistics is within our given time range, add it to
    % variable we'll use to calculate melt from a pulse. Times for this if
    % are given in myr.
    if statistics(i,2)/1e6 >= (maxyear+low) %&& statistics(i,2)/1e6 <= (maxyear+high)      
        meltf(i,1) = statistics(i,2);  %timestep

        %Subtract melt from previous timestep so we know how much melt is
        %produced in the given timestep, and then convert to m^3/ms. statistics 3
        %is the timestep size.
        meltf(i,2) = (statistics_smoothed(i,1) - statistics_smoothed(i-1,1))/(statistics(i,3)*year_to_second);

        % To compare if we had not smoothed the data
        meltff(i,1) = (statistics(i,16) - statistics(i-1,16))/(statistics(i,3)*year_to_second);%melt produced in m^3/ms
        
        %If we want to know what the original melt looked like.
        melt_orig(i,1) = statistics(i,16);
    end
end

%get rid of zeroes and shift time to zero.    
meltf = meltf(any(meltf,2),:);
meltff = meltff(any(meltff,2),:);

subplot(3,1,1)
scatter(meltf(:,1)/1e6,meltff(:,1),'.')
xlabel('Time (myr)');
ylabel('F (m^3/ms)');
title('Statistics melt flux')
set(gca,'FontSize',16)

subplot(3,1,2)
set(gca,'FontSize',16)
scatter(meltf(:,1)/1e6,meltff(:,1),'.')
hold on


%Find the minimum's between pulses, we'll then use this to determine where
%a pulse starts and ends, giving a time range to integrate for each pulse.
TF = islocalmin(meltf(:,2),'MinSeparation',2e2,'MinProminence', 8e-7);
meltmin(:,1) = meltf(TF,1);
meltmin(:,2) = meltf(TF,2);

plot(meltf(:,1)/1e6,meltf(:,2),'k','LineWidth',3)
hold on
scatter(meltmin(:,1)/1e6,meltmin(:,2),50,'ro','filled')


% Now we'll integrate each pulse by going from minimum to minimum. This can
% give us the total melt produced for each pulse.
for i=2:length(meltmin(:,1))

    % Get an equation for a line between the two minimums to isolate a pulse.
    [coefficients c] = polyfit([meltmin(i-1,1), meltmin(i,1)], [meltmin(i-1,2), meltmin(i,2)], 1);
    a = coefficients (1);
    b = coefficients (2);

    F = NaN;
    dt = NaN;
    t = NaN;
    yy = NaN;

    % Find time/F values and linear line based on our minimum line.
    for j=1:length(meltf(:,1))
        if(meltf(j,1) > meltmin(i-1,1) && meltf(j,1) < meltmin(i,1))
            % Melt after subtracting our line
            F = [F,meltf(j,2) - (meltf(j,1)*a+b)];
            % Timestep 
            dt = [dt, (meltf(j+1,1) - meltf(j-1,1))/2];
            % time
            t = [t, meltf(j,1)];
            % Our line for visualization
            yy = [yy, (meltf(j,1)*a+b)];
        elseif meltf(j,1) == meltmin(i,1)
            F = [F, meltf(j,2)- (meltf(j,1)*a+b)];
            dt = [dt, (meltf(j,1)- meltf(j-1,1))/2];
            t = [t, meltf(j,1)];
            yy = [yy, (meltf(j,1)*a+b)];
        elseif meltf(j,1) == meltmin(i-1,1)
            F = [F, meltf(j,2)- (meltf(j,1)*a+b)];
            dt = [dt, (meltf(j+1,1) - meltf(j,1))/2];
            t = [t, meltf(j,1)];
            yy = [yy, (meltf(j,1)*a+b)];
        end      
    end
    
    % Remove initial NaN
    F = rmmissing(F);
    yy = rmmissing(yy);
    dt = rmmissing(dt);
    t = rmmissing(t);
    
    %Plot area of each pulse
    subplot(3,1,2)
    set(gca,'FontSize',16)
    xlabel('Time (myr)');
    ylabel('F (m^3/ms)');
    title('Smoothed melt flux');
    plot(t/1e6,yy,'r','LineWidth',1.5);
    hold on

    subplot(3,1,3)
    set(gca,'FontSize',16)
    xlabel('Time (myr)');
    ylabel('F (m^3/ms)');
    title('Melt flux pulses');
    hold on

    % Shade every other pulse differently.
    if mod(i,2)==0
     area(t/1e6,F,'FaceAlpha',1,'EdgeAlpha', 0,'FaceColor',[0.6 0.8 1])
     plot(t/1e6,F,'k');
    else
     area(t/1e6,F,'FaceAlpha',1,'EdgeAlpha', 0,'FaceColor',[0.2 0.6 0.8])
     plot(t/1e6,F,'k');
    end
    
    % Multiply melt produced per time by dt to find total melt of the
    % pulse.
    for jj=1:length(dt)
        local_integration(jj) = F(jj)*dt(jj)*year_to_second;
    end  
    
   % Our averaged pulse size in km^3/km
   melt_integration(i-1) = sum(local_integration)*1e-6; % corrected to km^3/km
   flux_integration(i-1) = sum(local_integration)/((meltmin(i,1)-meltmin(i-1,1))*year_to_second);

   
   % And our average melt production in m^3/ms
   melt_average(i-1) = mean(F);
   melt_peak(i-1) = max(F);

   clear local_integration F t dt index yy
    
end

for i=1:length(melt_integration)
    if melt_integration(i) < 0
        melt_integration(i) = NaN;
        flux_integration(i) = NaN;
    end
end

melt_integration(isnan(melt_integration)) = [];
flux_integration(isnan(flux_integration)) = [];

pav = median(melt_integration);  
for ji = 1:length(melt_integration)
  melt_dif(ji) = pav - melt_integration(ji);
end

% Now put it all into a matrix for copying.
 final(1) = round(pav,1);
 final(3) = round(min(melt_integration),1);
 final(4) = round(max(melt_integration),1);
 %final(2)  = round(sqrt( mean( (melt_dif).^2 )),1 ); % For RMS on average

 final(5) = round(mean(flux_integration),7);
 final

clearvars -except final



