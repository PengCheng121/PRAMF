function [powratio2] = Rpower(f,p,orbit,sr,cutoff)

norbits = length(orbit);
nyquist=f(end);
nfpts = length(f);
powall = sum(p);
nfn_total = 0;

sr_num = length(sr);

for jj = 1:sr_num
    
    powratio = zeros(norbits,1);
    orbit_sed_p = 1./(orbit.*sr(jj)./100);
    orbit_max_min = zeros(2,norbits);
    
    ii = 0;
    
    for cf = 1:length(cutoff)
        nm = 0;
        ii = ii+1;
        
        if (orbit_sed_p(1)/cutoff(cf)) >= ((orbit_sed_p(2)-orbit_sed_p(1))/cutoff(cf))
            orbit_max_min(1,1) = orbit_sed_p(1) - (orbit_sed_p(2)-orbit_sed_p(1))/cutoff(cf);
            orbit_max_min(2,1) = orbit_sed_p(1) + (orbit_sed_p(2)-orbit_sed_p(1))/cutoff(cf);
        else
            orbit_max_min(1,1) = orbit_sed_p(1) - orbit_sed_p(1)/cutoff(cf);
            orbit_max_min(2,1) = orbit_sed_p(1) + orbit_sed_p(1)/cutoff(cf);
        end
        
        for i = 2:norbits-1
            left_cutoff = (orbit_sed_p(i)-orbit_sed_p(i-1))/cutoff(cf);
            right_cutoff = (orbit_sed_p(i+1)-orbit_sed_p(i))/cutoff(cf);
            if(left_cutoff > right_cutoff)
                orbit_max_min(1,i) = orbit_sed_p(i)-(orbit_sed_p(i+1)-orbit_sed_p(i))/cutoff(cf);
                orbit_max_min(2,i) = orbit_sed_p(i)+(orbit_sed_p(i+1)-orbit_sed_p(i))/cutoff(cf);
            else
                orbit_max_min(1,i) = orbit_sed_p(i)-(orbit_sed_p(i)-orbit_sed_p(i-1))/cutoff(cf);
                orbit_max_min(2,i) = orbit_sed_p(i)+(orbit_sed_p(i)-orbit_sed_p(i-1))/cutoff(cf);
            end
        end
        
        orbit_max_min(1,end) =  orbit_sed_p(end)-(orbit_sed_p(end)-orbit_sed_p(end-1))/cutoff(cf);
        orbit_max_min(2,end) =  orbit_sed_p(end)+(orbit_sed_p(end)-orbit_sed_p(end-1))/cutoff(cf);
        
        for i = 1 : norbits
            if nyquist <= orbit_sed_p(i) || f(2) >= orbit_sed_p(i)
                nm = nm + 1;
                continue;
            end
            
            nfmin=ceil(nfpts*orbit_max_min(1,i)/nyquist);
            
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
        
        powratio2(ii,1,jj) = cutoff(cf); % 截止频率大小
        powratio2(ii,2,jj) = sum(powratio)/powall;   % 天文旋回频段内能量比
        powratio2(ii,3,jj) = sum(powratio)/powall/sum(nfn_total);  % 天文旋回频段内能量比平均值
        powratio2(ii,4,jj) = (1 - sum(powratio)/powall);  % 噪声能量比
        powratio2(ii,5,jj) = (1 - sum(powratio)/powall)/(length(f) - sum(nfn_total));  % 噪声能量比平均值
        powratio2(ii,6,jj) = nm; % 无旋回的个数
    end
    
end

