function [h, p] = compare2Px2(group1, group2, group3, group4, conditions, ylbl, dim, markSz, color1, color2)
       
       if nargin < 7
           color1 = [0 0 0];
           color2 = [1 0 0];
       end
       if nargin < 6
           meanMarkSize = 15;
           markSize = 10;
       else
           markSize = markSz(1);
           meanMarkSize = markSz(2);
       end

       m = size(group1,1);
       m2 = size(group2,1);
       m3 = size(group3,1);
       m4 = size(group4,1);
       
       l_grey = [0.7 0.7 0.7];
       
       mean1 = nanmean(group1,1);
       std1 = sterr(group1,1);
       mean2 = nanmean(group2,1);
       std2 = sterr(group2,1);
       mean3 = nanmean(group3,1);
       std3 = sterr(group3,1);
       mean4 = nanmean(group4,1);
       std4 = sterr(group4,1);
       
       h = figure;
       line([1*ones(m,1)'; 2*ones(m2,1)'],[group1'; group2'], 'Color',l_grey); hold on;
       line([3*ones(m3,1)'; 4*ones(m4,1)'],[group3'; group4'], 'Color',l_grey); hold on;
       errorbar(1, mean1, std1,'LineStyle', 'none','LineWidth',1,'Color',color1,'CapSize',0,'Marker','.','MarkerSize',meanMarkSize);
       errorbar(2, mean2, std2,'LineStyle', 'none','LineWidth',1,'Color',color2,'CapSize',0,'Marker','.','MarkerSize',meanMarkSize);
       errorbar(3, mean3, std3,'LineStyle', 'none','LineWidth',1,'Color',color1,'CapSize',0,'Marker','.','MarkerSize',meanMarkSize);
       errorbar(4, mean4, std4,'LineStyle', 'none','LineWidth',1,'Color',color2,'CapSize',0,'Marker','.','MarkerSize',meanMarkSize);
       xlim([0.25 4.75]);
       ylim([0 inf])
       xticks([1:4]);
       xticklabels(conditions);
       ylabel(ylbl);
       [hp,p] = ttest(group1,group2);
       pt(5) = p;
       disp(p);
end

