% Script to find the viscosity ratio averages within running window
% First run runningwindow_periodicity on .mat from calculate_volumeflux, 
% then use periodicity and maxyear from that with this 

name = '_r2_tz31_v8_t550_p01.mat'
load(['../numerics/tilt/tilt', name])
load(['../numerics/viscosity_ratio/vr', name])

% Window size
interval = 25;

% alpha
signif = 0.01;

year_to_second = 365*24*60*60;


for j = 1:length(periodicity)
    
    % get all values within the time window.
    for i = 1:length(final410)
        if final410(i,1) >= (maxyear+(periodicity(j,2)-interval/2)) && final410(i,1) <= (maxyear+(periodicity(j,2)+interval/2))
            time200(i,1) = final410(i,1);
            time200(i,2) = final410(i,2);
            time200(i,3) = final410(i,4);
            time200(i,4) = final410(i,5);
            time200(i,5) = final410(i,7);
            time200(i,6) = final410(i,3);
            time200(i,7) = final410(i,11);
            time200(i,9) = final410(i,10);

            if length(final410(1,:)) == 11
                time200(i,8) = final410(i,11);
            else
                time200(i,8) = final410(i,12);
            end

        else
            time200(i,:) = 0;
        end
    end

    %Get rid of all the zero elements.
    time200 = time200(any(time200,2),:);

    if(isempty(time200) == 1)
        break
    end

    %Set first element as 0 so we can track over the total 205 Myr of data we
    %kept.
    time200(:,1) = time200(:,1)-time200(1,1);
   
%Integrate across variables   
for p=1:length(time200(:,1))
    if p==1
        dt(p,:) = (time200(p+1,1)-time200(p,1))/2;
    elseif p == length(time200)
        dt(p,:) = (time200(p,1) - time200(p-1,1))/2;
    else
        dt(p,:) = (time200(p+1,1) - time200(p-1,1))/2;
    end
end

%Sum up variables
    for jj=1:length(dt)
        vr_full(jj,:) = time200(jj,2)*dt(jj,:)*year_to_second;
        temp_full(jj,:) = time200(jj,4)*dt(jj,:)*year_to_second;
        vp_full(jj,:) = time200(jj,3)*dt(jj,:)*year_to_second;
        vm_full(jj,:) = time200(jj,6)*dt(jj,:)*year_to_second;
        atemp_full(jj,:) = time200(jj,7)*dt(jj,:)*year_to_second;
        ra_full(jj,:) = time200(jj,5)*dt(jj,:)*year_to_second;
        diff_full(jj,:) = time200(jj,8)*dt(jj,:)*year_to_second;
        mtemp_full(jj,:) = time200(jj,9)*dt(jj,:)*year_to_second;
    end
    

     periodicity(j,3) = sum(vr_full)/((time200(end,1)-time200(1,1))*year_to_second);
     periodicity(j,4) = sum(temp_full)/((time200(end,1)-time200(1,1))*year_to_second);
     periodicity(j,5) = sum(vp_full)/((time200(end,1)-time200(1,1))*year_to_second);
     periodicity(j,6) = sum(vm_full)/((time200(end,1)-time200(1,1))*year_to_second);
     periodicity(j,7) = sum(atemp_full)/((time200(end,1)-time200(1,1))*year_to_second);
     periodicity(j,8) = sum(ra_full)/((time200(end,1)-time200(1,1))*year_to_second);
     periodicity(j,10) = sum(diff_full)/((time200(end,1)-time200(1,1))*year_to_second);
     periodicity(j,11) = sum(mtemp_full)/((time200(end,1)-time200(1,1))*year_to_second);


    %%%%% Now do the same for tilt to get approximate pulse depth.
    for i = 1:length(tilt)
        if tilt(i,2) >= (maxyear+(periodicity(j,2)-interval/2)) && tilt(i,2) <= (maxyear+(periodicity(j,2)+interval/2))
            time300(i,1) = tilt(i,2);
            time300(i,2) = tilt(i,5);
        else
            time300(i,:) = 0;
        end
    end
    %Get rid of all the zero elements.
    time300 = time300(any(time300,2),:);

    if(isempty(time300) == 1)
        break
    end

    %Set first element as 0 so we can track over the total 205 Myr of data we
    %kept.
    time300(:,1) = time300(:,1)-time300(1,1);
   
%Integrate across variables   
for p=1:length(time300(:,1))
    if p==1
        dtt(p,:) = (time300(p+1,1)-time300(p,1))/2;
    elseif p == length(time300)
        dtt(p,:) = (time300(p,1) - time300(p-1,1))/2;
    else
        dtt(p,:) = (time300(p+1,1) - time300(p-1,1))/2;
    end
end

%Sum up variables
 for jj=1:length(dtt)
        depth_full(jj,:) = time300(jj,2)*dtt(jj,:)*year_to_second;
 end

     periodicity(j,9) = sum(depth_full)/((time300(end,1)-time300(1,1))*year_to_second);
     
     clearvars -except final410 periodicity maxyear interval year_to_second endtime tilt name2
     
end

for i = 1:length(periodicity(:,1))
    if (periodicity(i,1) == 0 || periodicity(i,3) == 0)
        periodicity(i,1) = 0;
        periodicity(i,2) = 0;
        periodicity(i,3) = 0;
        periodicity(i,4) = 0;
        periodicity(i,5) = 0;
        periodicity(i,6) = 0;
        periodicity(i,7) = 0;
        periodicity(i,8) = 0;
        periodicity(i,9) = 0;
        periodicity(i,10) = 0;
        periodicity(i,11) = 0;
    end
end

% Plot points
subplot(3,1,1)
scatter(periodicity(:,3),periodicity(:,1),200,'.')
title('Viscosity ratio vs. Periodicity');
xlabel('Viscosity ratio');
ylabel('Periodicity (Myr)');
grid on
hold on

subplot(3,1,2)
scatter(periodicity(:,7),periodicity(:,1),200,'.')
title('Plume non-adiabatic temperature vs. Periodicity');
xlabel('Nonadiabatic temperature (K)');
ylabel('Periodicity (Myr)');
grid on
hold on

subplot(3,1,3)
scatter(periodicity(:,2),periodicity(:,1),200,'.')
title('Time vs. Periodicity');
xlabel('Time after plume arrival (Myr)');
ylabel('Periodicity (Myr)');
grid on
hold on
