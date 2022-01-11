% First, load the volume flux .mat

cc1=0;

for yy = 50:1:50

interval = 25;
signif = 0.001;
%Find out when the max volume flux occurs
[M,I] = max(final410);

%yy = 200; %time to start from 105 for 100 Myr
maxsig = 50; %round((205-yy)/2,5);
xy = maxsig;  %xlimitdisp(['Most significant period is',num2str(1/f(jmax)),' with FAP of ',num2str(prob(jmax))])
%Set max year as the year max volume flux occurs
maxyear = final410(I(2),1)

%Make new matrix that only keeps data from 5 Myr before max volume flux,
%and 200 Myr after max volume flux. (Generally max is when plume first
%passes the chosen depth.)
for i = 1:length(final410)
    %if final410(i,1) >= (maxyear+(yy-interval/2)) && final410(i,1) <= (maxyear+(yy+interval/2))
    if final410(i,1) >= (maxyear + 70 - interval/2) && final410(i,1) <= (maxyear+200 + interval/2)
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
%plot(time200(:,1),time200(:,2))



%lombscargle(time200, 0);
 [P,f,alpha,M] = lomb(time200(:,2),time200(:,1));
 period = 1./f;
 
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
z = -log(1-(1-a).^(1/M));

for jj=1:length(periodp)
    if PyyP(jj) < z
        PyyP(jj) = NaN;
        periodp(jj,1) = NaN;
    end
end

periodp(any(isnan(periodp), 2), :) = [];
PyyP(any(isnan(PyyP), 2), :) = []; 
  
 mpp = periodp.*PyyP;
 mean2 = sum(mpp)/sum(PyyP)
 %mean3 = mean(sig_period)
% 
 
figure(1),
plot(time200(:,1),time200(:,2), 'k', 'LineWidth', 2)
title('Time vs. Volume flux 4 cm/yr (depth 100 km)');
xlabel('Time (Myr)');
ylabel('Volume flux (m^3/yr)');
xlim([0 400]);
% 
fig = figure
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 4.43]; 
plot(period,P, 'k', 'LineWidth', 1)
hold on
scatter(periodp,PyyP,'r*')
title('Lomb-normalized Periodogram')
ylabel('Power')
xlabel('Period (Myr)')
xlim([0 maxsig]);
xticks(0:10:maxsig)
%axes('Xlim', [0, 100], 'XTick', 0:5:100, 'NextPlot', 'add')
grid on;
    styles = {':','-.','--'};
    La = length(a);
    hold on;
    for i=1:La
        line([0,maxsig],[z(i),z(i)],'Color','k','LineStyle',styles{ceil(i*3/La)});
        text(maxsig-6,z(i)+0.03*max(P),strcat('\alpha = ',num2str(a(i))),'fontsize',10);
    end
%     



%ex = isempty('periodp');
if length(periodp) == 0
   peak_time(1,:) = 0;
else
   peak_time(1,:) = periodp;
end

format shortg
round(peak_time,2)
cc1 = [cc1,peak_time];
%round(peak_full,2)
%stop
%periodicity(yy,1) = peak_time;
%periodicity(yy,2) = yy;
%clearvars -except yy final410 periodicity peak_time
%clearvars -except cc1
end
% figure(3)
% plot(period,P,periodp,PyyP,'r*')
% xlim([0 xy]);

%clearvars -except periodicity
%clearvars -except peak_full peak_full5
%clearvars -except p peak_full1 peak_final peak_full ii jj info meannn maxxx minnn n44 n66
clear

