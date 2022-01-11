%%%Script to find the viscosity ratio averages over 25 Myr periods. From
%%%75
%%%Myr after plume arrival to 200 Myr%%%%
%%First run runningavg_periodicity on .mat from calculate_volumeflux, 
%% then use output from that and .mat output from calculate_RA... for this%%
%clear all
interval = 25;
signif = 0.001;
%Find out when the max volume flux occurs
[M,I] = max(final410);

%yy = 200; %time to start from 105 for 100 Myr
maxsig = 10; %round((205-yy)/2,5);
xy = maxsig;  %xlimitdisp(['Most significant period is',num2str(1/f(jmax)),' with FAP of ',num2str(prob(jmax))])
%Set max year as the year max volume flux occurs
maxyear = final410(I(2),1)
endtime = 200;
year_to_second = 365*24*60*60;



for j = 70:10:200 %time to start from 105 for 100 Myr
    
    %[M,I] = max(final410);
    maxsig = 12.5; %round((205-yy)/2,5);
    xy = maxsig;  %xlimitdisp(['Most significant period is',num2str(1/f(jmax)),' with FAP of ',num2str(prob(jmax))])
    %Set max year as the year max volume flux occurs
    %maxyear = final410(I(2),1);


    for i = 1:length(final410)
        if final410(i,1) >= (maxyear+(j-interval/2)) && final410(i,1) <= (maxyear+(j+interval/2))
            time200(i,1) = final410(i,1);
            time200(i,2) = final410(i,2);
            time200(i,3) = final410(i,4);
            time200(i,4) = final410(i,5);
            time200(i,5) = final410(i,7);
            time200(i,6) = final410(i,3);
            time200(i,7) = final410(i,9);
        else
            time200(i,:) = 0;
        end
    end
    %Get rid of all the zero elements.
    time200 = time200(any(time200,2),:);

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
    end
    

     periodicity(j,3) = sum(vr_full)/((time200(end,1)-time200(1,1))*year_to_second);
     periodicity(j,4) = sum(temp_full)/((time200(end,1)-time200(1,1))*year_to_second);
     periodicity(j,5) = sum(vp_full)/((time200(end,1)-time200(1,1))*year_to_second);
     periodicity(j,6) = sum(vm_full)/((time200(end,1)-time200(1,1))*year_to_second);
     periodicity(j,7) = sum(atemp_full)/((time200(end,1)-time200(1,1))*year_to_second);
     
     clearvars -except final410 periodicity maxyear interval year_to_second endtime
     
end

for i = 1:length(periodicity)
    if (periodicity(i,1) == 0 || periodicity(i,3) == 0)
        periodicity(i,1) = 0;
        periodicity(i,2) = 0;
        periodicity(i,3) = 0;
        periodicity(i,4) = 0;
        periodicity(i,5) = 0;
        periodicity(i,6) = 0;
        periodicity(i,7) = 0;
    end
    
    
end

periodicity = periodicity(any(periodicity,2),:);
save('per_70to200_25w_10s_41+200', 'periodicity')
clear

%round(ratio_average,2)
