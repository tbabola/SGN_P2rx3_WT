function makeDomHist(P)
lt_org = [255, 166 , 38]/255;
dk_org = [255, 120, 0]/255;
lt_blue = [50, 175, 242]/255;
dk_blue = [0, 13, 242]/255;

Lcounts = [];
Rcounts = [];
LpksAll = [];
RpksAll = [];
for j = 1 : length(P)
    Lpks = []; Rpks = [];
    for i = 1:length(P(j).locs)
        if P(j).L.mtrace(P(j).locs(i)) > P(j).R.mtrace(P(j).locs(i))
            Lpks = [Lpks; P(j).L.mtrace(P(j).locs(i))];
            LpksAll = [LpksAll; Lpks]
        else
            Rpks = [Rpks; P(j).R.mtrace(P(j).locs(i))];
            RpksAll = [RpksAll; Rpks]
        end
    end
    
    tMinutes = length(P(j).L.mtrace)/600;
    Lcounts = [Lcounts; histcounts(Lpks, [0:.07:.5])./tMinutes];
    Rcounts = [Rcounts; histcounts(Rpks, [0:.07:.5])./tMinutes];

end

Lmean = mean(Lcounts);
Rmean = mean(Rcounts);
Lsem = std(Lcounts)./sqrt(j);
Rsem = std(Rcounts)./sqrt(j);

hold on
L = barh(-Lmean);
L.FaceColor = lt_org; L.EdgeColor = lt_org;
R = barh(Rmean);
R.FaceColor = lt_blue; R.EdgeColor = lt_blue;
xlim([-2.5 2.5])
ylim([0 8])

for i = 1:7
    plot([-Lmean(i)-Lsem(i) -Lmean(i)+Lsem(i)],[i i], 'Color', dk_org, 'LineWidth', 2)
    plot([Rmean(i)-Rsem(i) Rmean(i)+Rsem(i)],[i i], 'Color', dk_blue, 'LineWidth', 2)
end
yticks([0 4 8])
yticklabels({0 25 50})
xticks([-2 -1 0 1 2])
xticklabels({2 1 0 1 2})

xlabel('Dominant events per minute')
ylabel('% \DeltaF/F')
fixAxis(gcf, gca)
end