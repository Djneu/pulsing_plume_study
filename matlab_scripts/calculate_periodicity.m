% First, load the volume flux .mat
clear
clf
cc1=0;

name = 'vf_r2_tz31_v8_t550_p1.mat'
load(['/home/djneuh/Desktop/movies/final_runs/data_files/volume_flux/newdata/', name])

% UMP initiation time in myr
umpt = 202;

% Alpha
signif = 0.01;

% Maximum period we would consider
maxsig = 50;
xy = maxsig; 

% Make new matrix that only keeps all data past UMP initiation.
for i = 2:length(final410)
    if final410(i,1) >= umpt
        time200(i,1) = final410(i,1);
        time200(i,2) = final410(i,2);
        timestep(i,1) = final410(i,1) - final410(i-1,1);
    else
        time200(i,:) = 0;
    end
end

%Get rid of all the zero elements.
time200 = time200(any(time200,2),:);

%Set first element as 0 so we can track over data we kept
time200(:,1) = time200(:,1)-time200(1,1);
time200(:,1) = smoothdata(time200(:,1),'movmean',2);

% Calculate periodogram
[P,f,alpha,MM] = lomb(time200(:,2),time200(:,1));

% Convert frequency to period.
period = 1./f;

% If period is outside our maximum range don't consider it.
for ii=1:length(period)
   if period(ii,1) >= maxsig
       P(ii,1) = NaN;
       period(ii,1) = NaN;
       alpha(ii,1) = NaN;
   end
end

period(any(isnan(period), 2), :) = [];
P(any(isnan(P), 2), :) = [];
alpha(any(isnan(alpha), 2), :) = [];
 
% Get rid of any signals that are not significant.
for iii=1:length(alpha)
    if alpha(iii,1) > signif
       sig_period(iii,1) = NaN;
       sig_P(iii,1) = NaN;
    else
       sig_period(iii,1) = period(iii);
       sig_P(iii,1) = P(iii);
    end
end

sig_period(any(isnan(sig_period), 2), :) = [];
sig_P(any(isnan(sig_P), 2), :) = []; 
 
Pyy = sig_P/max(sig_P);
mp = sig_period.*Pyy;
meann = sum(mp)/sum(Pyy);
  
  
TF = islocalmax(P,'MinSeparation',0); %'MinProminence', 5*10^8, 'MinSeparation',minutes(45),'SamplePoints',t  'MaxNumExtrema', 5
PyyP = P(TF);
periodp(:,1) = period(TF);

a = [signif];
z = -log(1-(1-a).^(1/MM));

for jj=1:length(periodp)
    if PyyP(jj) < z
        PyyP(jj) = NaN;
        periodp(jj,1) = NaN;
    end
end

periodp(any(isnan(periodp), 2), :) = [];
PyyP(any(isnan(PyyP), 2), :) = []; 
 
% Plot the volume flux range we considered.
subplot(2,1,1),
plot(time200(:,1),time200(:,2), 'k', 'LineWidth', 2)
title('Time vs. Volume flux');
xlabel('Time (Myr)');
ylabel('Volume flux (m^3/yr)');
xlim([0 400]);

% Plot the periodogram.
subplot(2,1,2)
plot(period,P, 'k', 'LineWidth', 1)
hold on
scatter(periodp,PyyP,'r*')
title('Lomb-normalized Periodogram')
ylabel('Power')
xlabel('Period (Myr)')
xlim([0 maxsig]);
xticks(0:10:maxsig)
grid on;
styles = {':','-.','--'};
La = length(a);
hold on;
for i=1:La
    line([0,maxsig],[z(i),z(i)],'Color','k','LineStyle',styles{ceil(i*3/La)});
    text(maxsig-6,z(i)+0.03*max(P),strcat('\alpha = ',num2str(a(i))),'fontsize',10);
end
     

if length(periodp) == 0
   peak_time(1,:) = 0;
else
   peak_time(1,:) = periodp;
end

format shortg
round(peak_time,1)
cc1 = [cc1,peak_time];

