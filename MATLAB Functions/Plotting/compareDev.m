function [h, p] = compareDev(group1, group2, group3, group4, conditions, ylbl, dim, markSz, color1, color2)
       
       if nargin < 9
           color1 = 'k';
           color2 = 'r';
       end
       if nargin < 7
           meanMarkSize = 15;
           markSize = 10;
       else
           markSize = markSz(1);
           meanMarkSize = markSz(2);
       end
       if size(group1,2) > size(group1,1)
           group1 = group1';
           group2 = group2';
           group3 = group3';
           group4 = group4';
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

       offset = 0.2;
       h = figure;
       xvalues = [1*ones(m,1); 2*ones(m2,1); 3*ones(m3,1);4*ones(m4,1)]-offset;
       xvalues = [xvalues [1*ones(m,1); 2*ones(m2,1); 3*ones(m3,1); 4*ones(m4,1)]+offset];
       plot(xvalues', [group1(:,1:2); group2(:,1:2); group3(:,1:2); group4(:,1:2)]','-','Color',l_grey);
       hold on;
       errorbar([1:4]-offset, [mean(group1(:,1)); mean(group2(:,1)); mean(group3(:,1)); mean(group4(:,1))], [sterr(group1(:,1),1); sterr(group2(:,1),1); sterr(group3(:,1),1); sterr(group4(:,1),1)],'LineStyle', 'none','LineWidth',1,'Color',color1,'CapSize',0,'Marker','.','MarkerSize',meanMarkSize);
       errorbar([1:4]'+offset, [mean(group1(:,2)); mean(group2(:,2)); mean(group3(:,2)); mean(group4(:,2))], [sterr(group1(:,2),1); sterr(group2(:,2),1); sterr(group3(:,2),1); sterr(group4(:,2),1)],'LineStyle', 'none','LineWidth',1,'Color',color2,'CapSize',0,'Marker','.','MarkerSize',meanMarkSize);
      
       xlim([0.5 4.5]);
       ylim([0 inf])
       xticks(1:4);
       xticklabels(conditions);
       ylabel(ylbl);
       handle = gcf;
       figQuality(gcf,gca,dim);
end

