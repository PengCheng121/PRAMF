clc;
clear;

close all;

tic;

%% load data
filename = "Newarkrank.csv";
load(filename)
Newarkrank(:,1) = Newarkrank(:,1)*0.3048;
% xx = 205:0.5:2000;
% xx = 2000:0.5:3600;
% xx = 205:0.5:3600;
xx = 2000:0.5:2750;
% xx = 2750:0.5:3600;
yy = interp1(Newarkrank(:,1),Newarkrank(:,2),xx,'linear');

dat = [xx',yy'];    % Two columns of data. The first colume is depth (Unit is m). The second column is value.

%% parameter setting
dt = dat(2,1)-dat(1,1);   % sampling rate of data (dat)
orbit = [405,133,101,42.8,34,25.6,21.5,20.2,17.5]; % target orbital parameters
% Newarkrank£¬202-204 Ma

red = 1;  % 0 = no remove red noise. else = robust AR1.
pad = 10000;   % zero-padding
srstep = 0.2; % step of sedimentation rates (Unit is cm/kyr).
sr1 = 0.2;   % begining sedimentation rates to be estimated (Unit is cm/kyr).
sr2 = 30;   % end sedimentation rates to be estimated (Unit is cm/kyr).
nsim = 600;   % number of Monte Carlo simulation
plotn = 1;  % 1 = plot power spectra. else = no plot
smoothwin = 0.25;   % smoothing window
linlog = 1;    % fit to S(f) or logS(f). 1 = linear; 2 = log
detrended = 1;    % 1 = detrending. else = no detrending
display = 1;    %  1 = show wait bar. else = no show.
method_powerspec = 1;    %  1 = Periodogram power spectral density estimate. 2 = Multitaper power spectral density estimate.
win = 100;    % Window length in detrending
method_smooth = 'lowess';    %  method of smoothing
method_cutoff = 1;  
  % 1 = Equally dividing method based on the number of intervals (The same number of intervals is used for all cycles in each Monte Carlo simulation). 
  % 2 = Equally dividing method based on the number of intervals (A different number of intervals is used for all cycles in each Monte Carlo simulation). 
  % 3 = Fixed frequency band method for the spatial domain spectrum besed on "findpeaks" function. In each Monte Carlo simulation, the frequency band is randomly generated between the maximum and minimum widths obtained by "findpeaks" function.
  % 4 = User-defined cutoff frequency limits (maximum fluctuation range of each period). 

number_interval_min = 2;  %  minimum number of intervals
number_interval_max = 20;   %  maximum number of intervals

if method_cutoff == 4
    cutoff_range = [40, 20, 8, 2, 2, 1, 0.5, 0.5, 1]; % The maximum fluctuation range of each period (Unit is kyr). For example, the two cutoff frequencies for the first cycle are calculated as follows: 1/(orbit(1)-cutoff_range(1)) and 1/(orbit(1)+cutoff_range(1)).
end

Unit_X_series = 'Depth (m)';
Unit_Y_series = 'Value';
Unit_X_spec = 'Frequency (cycle/m)';
Unit_Y_spec = 'Power';

%% detrending
if detrended == 1
    depth = dat(:,1);
    value = dat(:,2);
    data_detrend = smooth(depth,value,win/(depth(end)-depth(1)),method_smooth);
    dat2 = [depth,value-data_detrend];
else
    dat2 = dat;
end

%% power spectral analysis
if method_powerspec == 1
    [p,f] = periodogram(dat2(:,2),[],pad,1/dt);
else
    nw = 2;
    [p,f]=pmtm(dat2(:,2),nw,pad,1/dt);
end
P2 = p;
[pks,locs,w_pks,p_pks] = findpeaks(p,f);
w_pks_max = max(w_pks)/2;
w_pks_min = min(w_pks)/2;

%% remove red noise
[p,theored] = AR1noise(red,f,p,dt,smoothwin,linlog);
data = [f,p];

%% plot power spectra
if plotn == 1
    if red == 0
        figure;
        AX1 = subplot(2,1,1);
        plot(AX1,dat(:,1),dat(:,2),'b','LineWidth',1);
        hold on;
        if detrended == 1
            plot(AX1,dat2(:,1),dat2(:,2),'r','LineWidth',1);
            legend(AX1,'Data series','Detrended data series');
        else
            legend(AX1,'Data series');
        end
        title(AX1,'Data series');
        ylabel(AX1,Unit_Y_series);
        xlabel(AX1,Unit_X_series);
        set(AX1,'XMinorTick','on','YMinorTick','on');
        
        AX2 = subplot(2,1,2);
        plot(AX2,f,P2,'b','LineWidth',1);
        title(AX2,'Raw spectrum');
        xlabel(AX2,Unit_X_spec);
        ylabel(AX2,Unit_Y_spec);
        legend(AX2,'Power spectrum of data series');
        set(AX2,'XMinorTick','on','YMinorTick','on');
    else
        figure;
        AX1 = subplot(3,1,1);
        plot(AX1,dat(:,1),dat(:,2),'b','LineWidth',1);
        hold on;
        if detrended == 1
            plot(AX1,dat2(:,1),dat2(:,2),'r','LineWidth',1);
            legend(AX1,'Data series','Detrended data series');
        else
            legend(AX1,'Data series');
        end
        title(AX1,'Data series');
        ylabel(AX1,Unit_Y_series);
        xlabel(AX1,Unit_X_series);
        set(AX1,'XMinorTick','on','YMinorTick','on');
        
        AX2 = subplot(3,1,2);
        plot(AX2,f,P2,'b','LineWidth',1);
        hold on;
        plot(AX2,f,theored,'r','LineWidth',1);
        title(AX2,'Raw spectrum and red noise');
        xlabel(AX2,Unit_X_spec);
        ylabel(AX2,Unit_Y_spec);
        legend(AX2,'Power spectrum of data series','AR(1) noise');
        set(AX2,'XMinorTick','on','YMinorTick','on');
        
        AX3 = subplot(3,1,3);
        plot(AX3,f,p,'b','LineWidth',1);
        xlabel(AX3,Unit_X_spec);
        ylabel(AX3,Unit_Y_spec);
        title(AX3,'Red noise removed');
        legend(AX3,'Power spectrum of data series');
        set(AX3,'XMinorTick','on','YMinorTick','on');
    end
end

%% PRAMF analysis
sed_x = sr1:srstep:sr2;
mpts = length(sed_x);

nmi = zeros(mpts,nsim);
powratioi = zeros(mpts,nsim,4);
powratioi_minmax = zeros(mpts,nsim,4);
cutoff_nsim = zeros(nsim,1);

norbits = length(orbit);

if display == 1
    hwaitbar = waitbar(0,'Processing ... [CTRL + C to quit]','WindowStyle','normal','Name','Wait Bar');
end

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
        powall = sum(p);
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
            
            for i = 2:norbits
                orbit_max_min(1,i) = orbit_sed_p(i)-w_pks2;
                orbit_max_min(2,i) = orbit_sed_p(i)+w_pks2;
            end
            
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
    
    if display == 1
        waitbar(zz/nsim);
    end
    
end

if display == 1
    if ishandle(hwaitbar)
        close(hwaitbar);
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

%% plot
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

powyp = zeros(mpts,length(percent),3);

for i = 1:3
    powyp(:,:,i) = prctile(powratioi_minmax(:,:,i), percent,2);
    [nan_x,nan_y] = find(isnan(powyp(:,:,i)));
    powyp(nan_x,nan_y,i) = 0;
    nan_x = [];
    nan_y = [];
end

colorcode = [221/255,234/255,224/255; ...
    201/255,227/255,209/255; ...
    176/255,219/255,188/255;...
    126/255,201/255,146/255;...
    67/255,180/255,100/255];

sed_x = sed_x';

figure;
set(gcf,'unit','centimeters','position',[2,2,20,14])
set(gcf,'color','w');
for ii = 1:3
    ax1 = subplot(5,3,ii);
    for i = 1:npercent2
        fill(ax1,[sed_x; (fliplr(sed_x'))'],[powyp(:,npercent+1-i,ii); (fliplr(powyp(:,i,ii)'))'],colorcode(i,:),'LineStyle','none');
        hold on;
        leg{i} = num2str(percent(i));
    end
    plot(ax1,sed_x,powyp(:,npercent2+1,ii),'Color',[0,120/255,0],'LineWidth',1,'LineStyle','-');
    if ii ==3
        lgd = legend({'25%-75%','15.865%-84.135%','10%-90%','5%-95%','2.5%-97.5%'},'NumColumns',5,'FontSize',8,'Location','NorthEastOutside','FontName','Times New Roman','Position',[0.5 0.05 0 0]);
    end
    set(ax1,'XMinorTick','on','FontSize',8,'FontName','Times New Roman','Box','on');
    xlabel(ax1,'Sedimentation rate (cm/kyr)','FontSize',8,'FontName','Times New Roman');
    if ii == 1
        ylabel(ax1,{'Normalized', 'R_{Milankovitch}'},'FontSize',8,'FontName','Times New Roman');
        title('(a)','FontSize',8,'FontName','Times New Roman');
    elseif ii == 2
        ylabel(ax1,{'Normalized', 'R_{AveMi}'},'FontSize',8,'FontName','Times New Roman');
        title('(f)','FontSize',8,'FontName','Times New Roman');
    else
        ylabel(ax1,{'Normalized', '1-R_{AveNo}'},'FontSize',8,'FontName','Times New Roman');
        title('(k)','FontSize',8,'FontName','Times New Roman');
    end
    
    nmi2 = norbits - nmi;
    nmi2p = prctile(nmi2, percent,2);
    ax2 = subplot(5,3,ii+3);
    for i = 1:npercent2
        fill(ax2,[sed_x; (fliplr(sed_x'))'],[nmi2p(:,npercent+1-i); (fliplr(nmi2p(:,i)'))'],colorcode(i,:),'LineStyle','none');
        hold on;
    end
    plot(ax2,sed_x,nmi2p(:,npercent2+1),'Color',[0,120/255,0],'LineWidth',1,'LineStyle','-');
    ylabel(ax2,{'Number of', 'parameters'},'FontSize',8,'FontName','Times New Roman');
    set(ax2,'XMinorTick','on','FontSize',8,'FontName','Times New Roman','Box','on');
    xlabel(ax2,'Sedimentation rate (cm/kyr)','FontSize',8,'FontName','Times New Roman');
    if ii == 1
        title('(b)','FontSize',8,'FontName','Times New Roman');
    elseif ii == 2
        title('(g)','FontSize',8,'FontName','Times New Roman');
    else
        title('(l)','FontSize',8,'FontName','Times New Roman');
    end
    
    ax3 = subplot(5,3,ii+6);
    plot(ax3,sed_x,power_sig(:,ii),'r','LineWidth',1);
    if ii == 1
        ylabel(ax3,{'Accumulator of', 'R_{Milankovitch}'},'FontSize',8,'FontName','Times New Roman');
        title('(c)','FontSize',8,'FontName','Times New Roman');
    elseif ii == 2
        ylabel(ax3,{'Accumulator of', 'R_{AveMi}'},'FontSize',8,'FontName','Times New Roman');
        title('(h)','FontSize',8,'FontName','Times New Roman');
    else
        ylabel(ax3,{'Accumulator of', '1-R_{AveNo}'},'FontSize',8,'FontName','Times New Roman');
        title('(m)','FontSize',8,'FontName','Times New Roman');
    end
    set(ax3,'XMinorTick','on','FontSize',8,'FontName','Times New Roman','Box','on');
    xlabel(ax3,'Sedimentation rate (cm/kyr)','FontSize',8,'FontName','Times New Roman');
    
    ax4 = subplot(5,3,ii+9);
    plot(ax4,cutoff_nsim_sig(:,1,ii),cutoff_nsim_sig(:,2,ii),'ro','MarkerSize',3,'MarkerEdgeColor','b','MarkerFaceColor','b');
    xlim(ax4,[sed_x(1),sed_x(end)]);
    set(ax4,'XMinorTick','on','FontSize',8,'FontName','Times New Roman','Box','on');
    ylabel(ax4,{'Number of', 'intervals'},'FontSize',8,'FontName','Times New Roman');
    xlabel(ax4,'Sedimentation rate (cm/kyr)','FontSize',8,'FontName','Times New Roman');
    if ii == 1
        title('(d)','FontSize',8,'FontName','Times New Roman');
    elseif ii == 2
        title('(i)','FontSize',8,'FontName','Times New Roman');
    else
        title('(n)','FontSize',8,'FontName','Times New Roman');
    end
    
    ax5 = subplot(5,3,ii+12);
    plot(ax5,sed_x,powratioi_average(:,ii),'k','LineWidth',1);
    set(ax5,'XMinorTick','on','FontSize',8,'FontName','Times New Roman','Box','on');
    xlabel(ax5,'Sedimentation rate (cm/kyr)','FontSize',8,'FontName','Times New Roman');
    if ii == 1
        ylabel(ax5,{'Average','R_{Milankovitch}'},'FontSize',8,'FontName','Times New Roman');
        title('(e)','FontSize',8,'FontName','Times New Roman');
    elseif ii ==2
        ylabel(ax5,{'Average', 'R_{AveMi}'},'FontSize',8,'FontName','Times New Roman');
        title('(j)','FontSize',8,'FontName','Times New Roman');
    elseif ii ==3
        ylabel(ax5,{'Average', '1-R_{AveNo}'},'FontSize',8,'FontName','Times New Roman');
        title('(o)','FontSize',8,'FontName','Times New Roman');
    end
end


