clc;
clear;

close all;

tic;

%% load data
filename = "Newarkrank.csv";
load(filename)
Newarkrank(:,1) = Newarkrank(:,1)*0.3048;
xx = 205:0.5:2000;
% xx = 2000:0.5:3600;
% xx = 205:0.5:3600;
% xx = 2000:0.5:2750;
% xx = 2750:0.5:3600;
yy = interp1(Newarkrank(:,1),Newarkrank(:,2),xx,'linear');

dat = [xx',yy'];  % Two columns of data. The first colume is depth (Unit is m). The second column is value.

%% parameter setting
dt = dat(2,1)-dat(1,1);  % sampling rate of data (dat)
orbit = [405,133,101,42.8,34,25.6,21.5,20.2,17.5]; % Newarkrank£¬202-204 Ma

red = 1;  % 0 = no remove red noise. else = robust AR1.
pad = 10000;   % zero-padding
srstep = 0.5; % step of sedimentation rates (Unit is cm/kyr).
sr1 = 0.5;   % begining sedimentation rates to be estimated (Unit is cm/kyr).
sr2 = 30;   % end sedimentation rates to be estimated (Unit is cm/kyr).
nsim = 50;   % number of Monte Carlo simulation
smoothwin = 0.25;   % smoothing window
linlog = 1;    % fit to S(f) or logS(f). 1 = linear; 2 = log
detrended = 1;   % 1 = detrending. else = no detrending
display = 1;    %  1 = show wait bar. else = no show.
method_powerspec = 1;   %  1 = Periodogram power spectral density estimate. 2 = Multitaper power spectral density estimate.
win = 100;   % Window length in detrending
method_smooth = 'lowess';   %  method of smoothing
method_cutoff = 4;
  % 1 = Equally dividing method based on the number of intervals (The same number of intervals is used for all cycles in each Monte Carlo simulation). 
  % 2 = Equally dividing method based on the number of intervals (A different number of intervals is used for all cycles in each Monte Carlo simulation). 
  % 3 = Fixed frequency band method for the spatial domain spectrum besed on "findpeaks" function. In each Monte Carlo simulation, the frequency band is randomly generated between the maximum and minimum widths obtained by "findpeaks" function.
  % 4 = User-defined cutoff frequency limits (maximum fluctuation range of each period). 

window = 100;
step = 1;

number_interval_min = 2;  %  minimum number of intervals
number_interval_max = 20;   %  maximum number of intervals

if method_cutoff == 4
    cutoff_range = [40, 20, 8, 2, 2, 1, 0.5, 0.5, 1]; % The maximum fluctuation range of each period (Unit is kyr). For example, the two cutoff frequencies for the first cycle are calculated as follows: 1/(orbit(1)-cutoff_range(1)) and 1/(orbit(1)+cutoff_range(1)).
end

%% PRA analysis
% detrending
if detrended == 1
    depth = dat(:,1);
    value = dat(:,2);
    data_detrend = smooth(depth,value,win/(depth(end)-depth(1)),method_smooth);
    dat2 = [depth,value-data_detrend];
else
    dat2 = dat;
end

%% power spectral analysis
npts = fix(window/dt);
nrow = length(dat2(:,1));
m=fix((nrow-npts)/(step/dt));
sedrate = sr1:srstep:sr2;
sedrate = sedrate';
n = length(sedrate);

corrX = zeros(m,1);
epowratio_sig1 = zeros(m,n);
epowratio_sig2 = zeros(m,n);
epowratio_sig3 = zeros(m,n);
ecyclenum = zeros(m,n);

epowratio1 = zeros(m,n);
epowratio2 = zeros(m,n);
epowratio3 = zeros(m,n);

epowratioi_average1 = zeros(m,n);
epowratioi_average2 = zeros(m,n);
epowratioi_average3 = zeros(m,n);
norbits = length(orbit);

sed_x = sr1:srstep:sr2;
mpts = length(sed_x);


if display == 1
    hwaitbar = waitbar(0,'Processing ... [CTRL + C to quit]','WindowStyle','normal','Name','Wait Bar');
end

for iii = 1:m
    data = dat2(((iii-1)*(step/dt)+1):((iii-1)*(step/dt)+npts),:);
    
    
    if method_powerspec == 1
        [p,f] = periodogram(data(:,2),[],pad,1/dt);
    else
        nw = 2;
        [p,f]=pmtm(data(:,2),nw,pad,1/dt);
    end
    
    [pks,locs,w_pks,p_pks] = findpeaks(p,f);
    
    w_pks_max = max(w_pks)/2;
    w_pks_min = min(w_pks)/2;
    
    %% remove red noise
    [p,theored] = AR1noise(red,f,p,dt,smoothwin,linlog);
    % data = [f,p];
    
    %% PRA
    nmi = zeros(mpts,nsim);
    powratioi = zeros(mpts,nsim,4);
    powratioi_minmax = zeros(mpts,nsim,4);
    cutoff_nsim = zeros(nsim,1);
    powratio_var = zeros(mpts,nsim);
    
    powall = sum(p);
    
    for zz = 1:nsim
        
        j=1;
        
        if method_cutoff == 1
            cutoff = rand*(number_interval_max-number_interval_min)+number_interval_min;
            orbit_sed_p = 1./orbit;
            orbit_max_min = zeros(2,norbits);
            
            if (orbit_sed_p(1)/cutoff) >= ((orbit_sed_p(2)-orbit_sed_p(1))/cutoff)
                orbit_max_min(1,1) = orbit_sed_p(1) - (orbit_sed_p(2)-orbit_sed_p(1))/cutoff;
                orbit_max_min(2,1) = orbit_sed_p(1) + (orbit_sed_p(2)-orbit_sed_p(1))/cutoff;
            else
                orbit_max_min(1,1) = orbit_sed_p(1) - orbit_sed_p(1)/cutoff;
                orbit_max_min(2,1) = orbit_sed_p(1) + orbit_sed_p(1)/cutoff;
            end
            
            for i = 2:norbits-1
                left_cutoff = (orbit_sed_p(i)-orbit_sed_p(i-1))/cutoff;
                right_cutoff = (orbit_sed_p(i+1)-orbit_sed_p(i))/cutoff;
                if(left_cutoff > right_cutoff)
                    orbit_max_min(1,i) = orbit_sed_p(i)-(orbit_sed_p(i+1)-orbit_sed_p(i))/cutoff;
                    orbit_max_min(2,i) = orbit_sed_p(i)+(orbit_sed_p(i+1)-orbit_sed_p(i))/cutoff;
                else
                    orbit_max_min(1,i) = orbit_sed_p(i)-(orbit_sed_p(i)-orbit_sed_p(i-1))/cutoff;
                    orbit_max_min(2,i) = orbit_sed_p(i)+(orbit_sed_p(i)-orbit_sed_p(i-1))/cutoff;
                end
            end
            
            orbit_max_min(1,end) =  orbit_sed_p(end)-(orbit_sed_p(end)-orbit_sed_p(end-1))/cutoff;
            orbit_max_min(2,end) =  orbit_sed_p(end)+(orbit_sed_p(end)-orbit_sed_p(end-1))/cutoff;
        end
        
        if method_cutoff == 2
            
            for i = 1:norbits
                cutoff(i) = rand*(number_interval_max-number_interval_min)+number_interval_min;
            end
            
            orbit_sed_p = 1./orbit;
            orbit_max_min = zeros(2,norbits);
            
            if (orbit_sed_p(1)/cutoff(1)) >= ((orbit_sed_p(2)-orbit_sed_p(1))/cutoff(1))
                orbit_max_min(1,1) = orbit_sed_p(1) - (orbit_sed_p(2)-orbit_sed_p(1))/cutoff(1);
                orbit_max_min(2,1) = orbit_sed_p(1) + (orbit_sed_p(2)-orbit_sed_p(1))/cutoff(1);
            else
                orbit_max_min(1,1) = orbit_sed_p(1) - orbit_sed_p(1)/cutoff(1);
                orbit_max_min(2,1) = orbit_sed_p(1) + orbit_sed_p(1)/cutoff(1);
            end
            
            for i = 2:norbits-1
                left_cutoff = (orbit_sed_p(i)-orbit_sed_p(i-1))/cutoff(i);
                right_cutoff = (orbit_sed_p(i+1)-orbit_sed_p(i))/cutoff(i);
                if(left_cutoff > right_cutoff)
                    orbit_max_min(1,i) = orbit_sed_p(i)-(orbit_sed_p(i+1)-orbit_sed_p(i))/cutoff(i);
                    orbit_max_min(2,i) = orbit_sed_p(i)+(orbit_sed_p(i+1)-orbit_sed_p(i))/cutoff(i);
                else
                    orbit_max_min(1,i) = orbit_sed_p(i)-(orbit_sed_p(i)-orbit_sed_p(i-1))/cutoff(i);
                    orbit_max_min(2,i) = orbit_sed_p(i)+(orbit_sed_p(i)-orbit_sed_p(i-1))/cutoff(i);
                end
            end
            
            orbit_max_min(1,end) =  orbit_sed_p(end)-(orbit_sed_p(end)-orbit_sed_p(end-1))/cutoff(end);
            orbit_max_min(2,end) =  orbit_sed_p(end)+(orbit_sed_p(end)-orbit_sed_p(end-1))/cutoff(end);
        end
        
        if method_cutoff == 3
            w_pks_rand = (rand*(w_pks_max-w_pks_min)+w_pks_min);
        end
        
        if method_cutoff == 4
            
            for i = 1:norbits
                
                cutoff = rand;
                
                orbit_max_min(1,i) = 1./(orbit(i)+cutoff_range(i)*cutoff);
                orbit_max_min(2,i) = 1./(orbit(i)-cutoff_range(i)*cutoff);
                orbit_sed_p = 1./orbit;
            end
            
        end

        for sr = sr1:srstep:sr2
            
            nm = 0;
            
            powratio = zeros(norbits,1);

            nfn_total = 0;
            
            F = f.*(sr./100);
            nyquist=F(end);
            nfpts = length(F);
            
            if method_cutoff == 3
                
                w_pks2 = w_pks_rand.*(sr./100);
                
                orbit_sed_p = 1./orbit;
                orbit_max_min = zeros(2,norbits);
                
                if orbit_sed_p(1) >= w_pks2
                    orbit_max_min(1,1) = orbit_sed_p(1) - w_pks2;
                    orbit_max_min(2,1) = orbit_sed_p(1) + w_pks2;
                else
                    orbit_max_min(1,1) = 0;
                    orbit_max_min(2,1) = orbit_sed_p(1) + orbit_sed_p(1);
                end
                
                for i = 2:norbits-1
                    orbit_max_min(1,i) = orbit_sed_p(i)-w_pks2;
                    orbit_max_min(2,i) = orbit_sed_p(i)+w_pks2;
                end
                
                orbit_max_min(1,end) =  orbit_sed_p(end)-w_pks2;
                orbit_max_min(2,end) =  orbit_sed_p(end)+w_pks2;
            end
            
            indicate_value = zeros(1,norbits-1);
            for i = 1:norbits-1
                if (orbit_max_min(2,i) - orbit_max_min(1,i+1)) >= 0
                    indicate_value(1,i) = 1;
                end
            end
            
            
            if sum(indicate_value) > 0
                nm = norbits;
            else
                for i = 1 : norbits
                    
                    if nyquist <= orbit_sed_p(i) || F(2) >= orbit_sed_p(i)
                        nm = nm + 1;
                        continue;
                    end
                    
                    nfmin=ceil(nfpts*orbit_max_min(1,i)/nyquist);
                    
                    if nfmin == 0
                        nfmin = 1;
                    end
                    
                    nfmax=fix(nfpts*orbit_max_min(2,i)/nyquist);
                    
                    if orbit_max_min(2,i)/nyquist>1
                        nfmax = nfpts;
                    end
                    
                    nfn=nfmax-nfmin+1;
                    spq=zeros(1,nfn);
                    ij=1;
                    for q=nfmin:nfmax
                        spq(ij)=p(q);
                        ij=ij+1;
                    end
                    
                    powratio(i,1) = sum(spq);
                    
                    if sum(spq) <= 0
                        nm = nm + 1;
                    end
                    
                    nfn_total(i,1) = nfn;
                    
                end
            end
            
            powratioi(j,zz,1) = sum(powratio)/powall;
            powratioi(j,zz,2) = sum(powratio)/powall/sum(nfn_total);
            powratioi(j,zz,3) = 1 - (1 - sum(powratio)/powall)/(length(f) - sum(nfn_total));
            
            nmi(j,zz) = nm;
            
            powratio_max = max(powratio);
            for i = 1 : norbits
                powratio_var(j,zz) = powratio_var(j,zz) + abs(powratio(i,1) - powratio_max);
            end
            powratio_var(j,zz) = powratio_var(j,zz)/powratioi(j,zz,1);
            
            j=j+1;
        end
        
        for i = 1:3
            powratioi_minmax(:,zz,i) = (powratioi(:,zz,i)-min(powratioi(:,zz,i))) ./ (max(powratioi(:,zz,i))-min(powratioi(:,zz,i)));
        end
        
        if method_cutoff == 1 || method_cutoff == 4
            cutoff_nsim(zz,1) = cutoff;
        end
        
        if method_cutoff == 2
            cutoff_nsim(zz,1) = mean(cutoff);
        end
        
        if method_cutoff == 3
            cutoff_nsim(zz,1) = w_pks2;
        end
    end
    
    power_sig = zeros(mpts,3);
    cutoff_nsim_sig = zeros(nsim,2,3);
    powratioi_average = zeros(mpts,3);
    
    for i = 1:3
        
        powratioi_average(:,i) = mean(powratioi(:,:,i),2);
        
        for zz = 1:nsim
            if max(powratioi_minmax(:,zz,i)) ==1
                index = find(powratioi_minmax(:,zz,i) == 1);
                for jj = 1:length(index)
                    power_sig(index(jj),i) = power_sig(index(jj),i)+1;
                    cutoff_nsim_sig(zz,1,i) = sed_x(1,index(jj));
                    cutoff_nsim_sig(zz,2,i) = cutoff_nsim(zz,1);
                end
            end
        end
    end
    power_sig = power_sig/nsim;
    
    checkbox1median = 50;
    checkbox50 = [25 75];
    checkbox68 = [15.865 84.135];
    checkbox80 = [10 90];
    checkbox90 = [5 95];
    checkbox95 = [2.5 97.5];
    
    percent = [checkbox1median, checkbox50, checkbox68, ...
        checkbox80,checkbox90, checkbox95];

    percent = sort(unique(percent));
    percent = percent(percent~=0);
    if isempty(percent)
        percent = 50;
    end
    npercent  = length(percent);
    npercent2 = (length(percent)-1)/2;
    
    for i = 1:3
        powyp(:,:,i) = prctile(powratioi_minmax(:,:,i), percent,2);
    end
    nmi2 = norbits - nmi;
    nmi2p = prctile(nmi2, percent,2);
    
    corrX(iii,1) = dat2(fix(npts/2+(iii-1)*(step/dt)),1);
    
    epowratio1(iii,:) = powyp(:,npercent2+1,1);
    epowratio2(iii,:) = powyp(:,npercent2+1,2);
    epowratio3(iii,:) = powyp(:,npercent2+1,3);
    
    epowratio_sig1(iii,:) = power_sig(:,1);
    epowratio_sig2(iii,:) = power_sig(:,2);
    epowratio_sig3(iii,:) = power_sig(:,3);
    ecyclenum(iii,:) = nmi2p(:,npercent2+1);
    
    epowratioi_average1(iii,:) = powratioi_average(:,1);
    epowratioi_average2(iii,:) = powratioi_average(:,2);
    epowratioi_average3(iii,:) = powratioi_average(:,3);
    
    if display == 1
        waitbar(iii/m);
    end
end


if display == 1
    if ishandle(hwaitbar)
        close(hwaitbar);
    end
end

for i = 1:m
    corrX(i,1) = dat2(fix(npts/2+(i-1)*(step/dt)),1);
end

figure
set(gcf,'unit','centimeters','position',[10,5,18,10])
set(gcf,'color','w');
ax1 = subplot(1,4,1);
[C,h]=contour(sedrate,corrX,epowratio1);
h.Fill = 'on';
u1 = colorbar('southoutside','FontSize',8,'FontName','Times New Roman');
u1.Label.String = 'Median R_{Milankovitch}';
colormap(jet)
shading interp
ylabel('Depth (m)','FontSize',8,'FontName','Times New Roman')
xlabel('Sedimentation rate (cm/kyr)','FontSize',8,'FontName','Times New Roman')
zlabel('powratioi','FontSize',8,'FontName','Times New Roman')
figurename ='Median R_{Milankovitch}';
title(figurename)
set(ax1,'Ydir','reverse','FontSize',8,'FontName','Times New Roman');

ax2 = subplot(1,4,2);
[C,h]=contour(sedrate,corrX,epowratio2);
h.Fill = 'on';
u1 = colorbar('southoutside','FontSize',8,'FontName','Times New Roman');
u1.Label.String = 'Median R_{AveMi}';
colormap(jet)
shading interp
ylabel('Depth (m)','FontSize',8,'FontName','Times New Roman')
xlabel('Sedimentation rate (cm/kyr)','FontSize',8,'FontName','Times New Roman')
zlabel('powratioi','FontSize',8,'FontName','Times New Roman')
figurename ='Median R_{AveMi}';
title(figurename)
set(ax2,'Ydir','reverse','FontSize',8,'FontName','Times New Roman');

ax3 = subplot(1,4,3);
[C,h]=contour(sedrate,corrX,epowratio3);
h.Fill = 'on';
u1 = colorbar('southoutside','FontSize',8,'FontName','Times New Roman');
u1.Label.String = 'Median 1-R_{AveNo}';
colormap(jet)
shading interp
ylabel('Depth (m)','FontSize',8,'FontName','Times New Roman')
xlabel('Sedimentation rate (cm/kyr)','FontSize',8,'FontName','Times New Roman')
zlabel('powratioi','FontSize',8,'FontName','Times New Roman')
figurename ='Median 1-R_{AveNo}';
title(figurename)
set(ax3,'Ydir','reverse','FontSize',8,'FontName','Times New Roman');

ax4 = subplot(1,4,4);
zlevs = 0:1:length(orbit);
[C,h]=contour(sedrate,corrX,ecyclenum,zlevs);
h.Fill = 'on';
colormap(jet)
shading interp
u3 = colorbar('southoutside','FontSize',8,'FontName','Times New Roman');
u3.Label.String = 'Number of parameters';
xlabel('Sedimentation rate (cm/kyr)','FontSize',8,'FontName','Times New Roman')
zlabel('#')
figurename = 'Number of parameters';
title(figurename)
set(ax4,'Ydir','reverse','FontSize',8,'FontName','Times New Roman');


figure
set(gcf,'unit','centimeters','position',[10,5,18,10])
set(gcf,'color','w');
ax1 = subplot(1,4,1);
[C,h]=contour(sedrate,corrX,epowratio_sig1);
h.Fill = 'on';
u1 = colorbar('southoutside','FontSize',8,'FontName','Times New Roman');
u1.Label.String = 'Accumulator of R_{Milankovitch}';
colormap(jet)
shading interp
ylabel('Depth (m)','FontSize',8,'FontName','Times New Roman')
xlabel('Sedimentation rate (cm/kyr)','FontSize',8,'FontName','Times New Roman')
zlabel('powratioi','FontSize',8,'FontName','Times New Roman')
figurename ='Accumulator of R_{Milankovitch}';
title(figurename)
set(ax1,'Ydir','reverse','FontSize',8,'FontName','Times New Roman');

ax2 = subplot(1,4,2);
[C,h]=contour(sedrate,corrX,epowratio_sig2);
h.Fill = 'on';
u1 = colorbar('southoutside','FontSize',8,'FontName','Times New Roman');
u1.Label.String = 'Accumulator of R_{AveMi}';
colormap(jet)
shading interp
ylabel('Depth (m)','FontSize',8,'FontName','Times New Roman')
xlabel('Sedimentation rate (cm/kyr)','FontSize',8,'FontName','Times New Roman')
zlabel('powratioi','FontSize',8,'FontName','Times New Roman')
figurename ='Accumulator of R_{AveMi}';
title(figurename)
set(ax2,'Ydir','reverse','FontSize',8,'FontName','Times New Roman');

ax3 = subplot(1,4,3);
[C,h]=contour(sedrate,corrX,epowratio_sig3);
h.Fill = 'on';
u1 = colorbar('southoutside','FontSize',8,'FontName','Times New Roman');
u1.Label.String = 'Accumulator of 1-R_{AveNo}';
colormap(jet)
shading interp
ylabel('Depth (m)','FontSize',8,'FontName','Times New Roman')
xlabel('Sedimentation rate (cm/kyr)','FontSize',8,'FontName','Times New Roman')
zlabel('powratioi','FontSize',8,'FontName','Times New Roman')
figurename ='Accumulator of 1-R_{AveNo}';
title(figurename)
set(ax3,'Ydir','reverse','FontSize',8,'FontName','Times New Roman');

ax4 = subplot(1,4,4);
zlevs = 0:1:length(orbit);
[C,h]=contour(sedrate,corrX,ecyclenum,zlevs);
h.Fill = 'on';
colormap(jet)
shading interp
u3 = colorbar('southoutside','FontSize',8,'FontName','Times New Roman');
u3.Label.String = 'Number of parameters';
xlabel('Sedimentation rate (cm/kyr)','FontSize',8,'FontName','Times New Roman')
zlabel('#')
figurename = 'Number of parameters';
title(figurename)
set(ax4,'Ydir','reverse','FontSize',8,'FontName','Times New Roman');

figure
set(gcf,'unit','centimeters','position',[10,5,18,10])
set(gcf,'color','w');
ax1 = subplot(1,4,1);
[C,h]=contour(sedrate,corrX,epowratioi_average1);
h.Fill = 'on';
u1 = colorbar('southoutside','FontSize',8,'FontName','Times New Roman');
u1.Label.String = 'Average R_{Milankovitch}';
colormap(jet)
shading interp
ylabel('Depth (m)','FontSize',8,'FontName','Times New Roman')
xlabel('Sedimentation rate (cm/kyr)','FontSize',8,'FontName','Times New Roman')
zlabel('powratioi','FontSize',8,'FontName','Times New Roman')
figurename ='Average R_{Milankovitch}';
title(figurename)
set(ax1,'Ydir','reverse','FontSize',8,'FontName','Times New Roman');

ax2 = subplot(1,4,2);
[C,h]=contour(sedrate,corrX,epowratioi_average2);
h.Fill = 'on';
u1 = colorbar('southoutside','FontSize',8,'FontName','Times New Roman');
u1.Label.String = 'Average R_{AveMi}';
colormap(jet)
shading interp
ylabel('Depth (m)','FontSize',8,'FontName','Times New Roman')
xlabel('Sedimentation rate (cm/kyr)','FontSize',8,'FontName','Times New Roman')
zlabel('powratioi','FontSize',8,'FontName','Times New Roman')
figurename ='Average R_{AveMi}';
title(figurename)
set(ax2,'Ydir','reverse','FontSize',8,'FontName','Times New Roman');

ax3 = subplot(1,4,3);
[C,h]=contour(sedrate,corrX,epowratioi_average3);
h.Fill = 'on';
u1 = colorbar('southoutside','FontSize',8,'FontName','Times New Roman');
u1.Label.String = 'Average 1-R_{AveNo}';
colormap(jet)
shading interp
ylabel('Depth (m)','FontSize',8,'FontName','Times New Roman')
xlabel('Sedimentation rate (cm/kyr)','FontSize',8,'FontName','Times New Roman')
zlabel('powratioi','FontSize',8,'FontName','Times New Roman')
figurename ='Average 1-R_{AveNo}';
title(figurename)
set(ax3,'Ydir','reverse','FontSize',8,'FontName','Times New Roman');

ax4 = subplot(1,4,4);
zlevs = 0:1:length(orbit);
[C,h]=contour(sedrate,corrX,ecyclenum,zlevs);
h.Fill = 'on';
colormap(jet)
shading interp
u3 = colorbar('southoutside','FontSize',8,'FontName','Times New Roman');
u3.Label.String = 'Number of parameters';
xlabel('Sedimentation rate (cm/kyr)','FontSize',8,'FontName','Times New Roman')
zlabel('#')
figurename = 'Number of parameters';
title(figurename)
set(ax4,'Ydir','reverse','FontSize',8,'FontName','Times New Roman');

