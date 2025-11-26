function [nm,powratio2,powratio] = ratiovalue(f,data_power,orbit,sr,method_cutoff,cutoff,cutoff_value)

nm = 0;
norbits = length(orbit);

powratio = zeros(norbits,1);
powall = sum(data_power);
nfn_total = 0;

F = f.*(sr./100);
nyquist=F(end);
nfpts = length(F);



if method_cutoff == 1

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
    
    
elseif method_cutoff == 2
    
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
    
elseif method_cutoff == 3
    
%     [~,~,w_pks,~] = findpeaks(data_power,F);
%     w_pks = mean(w_pks)/2;
    
    w_pks = cutoff_value.*(sr./100);
%     w_pks = cutoff_value;
        
    orbit_sed_p = 1./orbit;
    orbit_max_min = zeros(2,norbits);
    
    
    if orbit_sed_p(1) >= w_pks
        orbit_max_min(1,1) = orbit_sed_p(1) - w_pks;
        orbit_max_min(2,1) = orbit_sed_p(1) + w_pks;
    else
        orbit_max_min(1,1) = 0;
        orbit_max_min(2,1) = orbit_sed_p(1) + orbit_sed_p(1);
    end
    
    
    
    for i = 2:norbits-1
        orbit_max_min(1,i) = orbit_sed_p(i)-w_pks;
        orbit_max_min(2,i) = orbit_sed_p(i)+w_pks;
    end
    
    orbit_max_min(1,end) =  orbit_sed_p(end)-w_pks;
    orbit_max_min(2,end) =  orbit_sed_p(end)+w_pks;
    
    

    
else

    orbit_max_min = 1./(cutoff);
    orbit_sed_p = 1./orbit;

end


% figure;
% set(gcf,'unit','centimeters','position',[10,5,7.5,7.5])
% subplot('Position',[0.15 0.7 0.55 0.25]);
% plot(F,data_power);
% hold on;
% for i = 1:norbits
%     xline(orbit_max_min(1,i));
%     xline(orbit_max_min(2,i));
%     hold on;
% end
% xlim([0,0.08]);
% ylabel('Power','FontSize',8,'FontName','Times New Roman');
% xlabel('Frequency (cycles/m)','FontSize',8,'FontName','Times New Roman');


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
            spq(ij)=data_power(q);
            ij=ij+1;
        end
        
        powratio(i,1) = sum(spq)/powall;
        
        if sum(spq) <= 0
            nm = nm + 1;
        end
        
        nfn_total(i,1) = nfn;
        
    end
    
end

powratio2 = sum(powratio);
% powratio2 = (1 - sum(powratio)/powall)/(length(f) - sum(nfn_total)); 









