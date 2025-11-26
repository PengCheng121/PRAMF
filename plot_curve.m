function plot_curve(powratio2,cutoff,sr,sr_optimal,parameter)

plot_noise = zeros(1,length(cutoff));

figure;
for ii = 1:length(sr)
    for jj = 1:length(cutoff)
        plot_noise(ii,jj) = powratio2(jj,parameter,ii);
    end
    if(sr(ii) == sr_optimal)
        plot(cutoff,plot_noise(ii,:),'r','LineWidth',1);
        hold on;
    else
        if (sr(ii) < sr_optimal)
            plot(cutoff,plot_noise(ii,:),'b','LineStyle','--');
            hold on;
        else
            plot(cutoff,plot_noise(ii,:),'k','LineStyle','--');
            hold on;
        end
    end
end