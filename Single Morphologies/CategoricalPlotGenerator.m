% CategoricalPlotGenerator.m
% Description: Plot categorical scatter plot of distributions of different morphological
%              regions (fully automatic - change title labels if desired)
% Usage: Run using MATLAB in the same folder as the LFPy script output
% Author(s): Ilhan Bok
% Last Modified: Jan. 19, 2022

close all; clear; clc;

parseRepr2;

% Whether to use absolute magnetization value (vs. signed averages) in
% comparing neuron types (recommended: 1)
compareABSVAL = 1;

% Correction factor to scale cell magnetization by
sf = 1;

% CHANGE THESE IF DESIRED
groups = {
'Middle Temporal Gyrus Layer 3'
'Middle Temporal Gyrus Layer 6'
'Frontal Lobe Layer 3'
'Middle Frontal Gyrus 3'
};

types = {
    'MTG3'
    'MTG6'
    'FLL3'
    'MFG3'
};
means = [];
for g = (0:3)
    close all;
    for i = (1:4)
        eval(['soma' num2str(i) '=cell' num2str(4*g + i) '.soma''']);
    end
    % Assemble soma data
    gs1 = cell(1,length(soma1));
    gs1(:) = {'Soma (Cell I)'};

    gs2 = cell(1,length(soma2));
    gs2(:) = {'Soma (Cell II)'};

    gs3 = cell(1,length(soma3));
    gs3(:) = {'Soma (Cell III)'};

    gs4 = cell(1,length(soma4));
    gs4(:) = {'Soma (Cell IV)'};

    somag = [soma1./length(soma1) soma2./length(soma2) soma3./length(soma3) ...
        soma4./length(soma4)];
    somag = somag.*length(somag)./4;
    somag = [logmean(soma1) logmean(soma2) logmean(soma3) logmean(soma4)];

    gsm = cell(1,length(somag));
    gsm(:) = {'Soma (Mean)'};

    for i = (1:4)
        try
            eval(['axon' num2str(i) '=cell' num2str(4*g + i) '.axon''']);
        catch Exception
            disp('Cell has no axon (it''s fine');
            eval(['axon' num2str(i) '=[]']);
        end
    end
    % Assemble axon data
    ga1 = cell(1,length(axon1));
    ga1(:) = {'Axon (Cell I)'};

    ga2 = cell(1,length(axon2));
    ga2(:) = {'Axon (Cell II)'};

    ga3 = cell(1,length(axon3));
    ga3(:) = {'Axon (Cell III)'};

    ga4 = cell(1,length(axon4));
    ga4(:) = {'Axon (Cell IV)'};

    axonag = [axon1./length(axon1) axon2./length(axon2) axon3./length(axon3) ...
        axon4./length(axon4)];
    axonag = axonag.*length(axonag)./4;
    
    axonag = [];
    for i = (1:4)
        if ~isempty(eval(['axon' num2str(i)]))
            axonag(end+1) = eval(['logmean(axon' num2str(i) ')'])
        end
    end

    gam = cell(1,length(axonag));
    gam(:) = {'Axon (Mean)'};

    for i = (1:4)
        eval(['dend' num2str(i) '=cell' num2str(4*g + i) '.dend''']);
    end
    
    gd1 = cell(1,length(dend1));
    gd1(:) = {'Dendrites (Cell I)'};

    gd2 = cell(1,length(dend2));
    gd2(:) = {'Dendrites (Cell II)'};

    gd3 = cell(1,length(dend3));
    gd3(:) = {'Dendrites (Cell III)'};

    gd4 = cell(1,length(dend4));
    gd4(:) = {'Dendrites (Cell IV)'};

    dendag = [dend1./length(dend1) dend2./length(dend2) dend3./length(dend3) ...
        dend4./length(dend4)];
    dendag = dendag.*length(dendag)./4;
    dendag = [logmean(dend1) logmean(dend2) logmean(dend3) logmean(dend4)];

    gdm = cell(1,length(dendag));
    gdm(:) = {'Dendrites (Mean)'};
    
    allag = [soma1 soma2 soma3 soma4 axon1 axon2 axon3 axon4 dend1 dend2 dend3 dend4];
    means(end+1) = mean(abs([soma1 axon1 dend1]));
    means(end+1) = mean(abs([soma2 axon2 dend2]));
    means(end+1) = mean(abs([soma3 axon3 dend3]));
    means(end+1) = mean(abs([soma4 axon4 dend4]));
    
    gtm = cell(1,length(allag));
    gtm(:) = {'Cell (Aggregate)'};

    merge = {gtm, gsm, gs1, gs2, gs3, gs4, gam, ga1, ga2, ga3, ga4, gdm, gd1, gd2, gd3, gd4};

    gr = horzcat(merge{:});
    f = 0.7;
     if compareABSVAL
        x = cat(2,abs(allag),abs(somag),abs(soma1),abs(soma2),abs(soma3),abs(soma4),...
        abs(axonag),abs(axon1),abs(axon2),abs(axon3),abs(axon4),...
        abs(dendag),abs(dend1),abs(dend2),abs(dend3),abs(dend4));
     else
        x = cat(2,allag,somag,soma1,soma2,soma3,soma4,...
        axonag,axon1,axon2,axon3,axon4,...
        dendag,dend1,dend2,dend3,dend4);
     end
    eval(['allagt' num2str(g+1) '=allag']);
    eval(['somagt' num2str(g+1) '=somag']);
    eval(['axonagt' num2str(g+1) '=axonag']);
    eval(['dendagt' num2str(g+1) '=dendag']);

    f = 0.7;
    fig = figure;
    switch(types{g+1})
        case 'MTG3'
            %c = 'kgggggbbbbrrrrr';
            c = [f f f;
                 f/2 1 f/2;
                 f 1 f;
                 f 1 f;
                 f 1 f;
                 f 1 f;
                 f/2 f/2 1;
                 f f 1;
                 f f 1;
                 f f 1;
                 1 f/2 f/2;
                 1 f f;
                 1 f f;
                 1 f f
                 1 f f];
        case 'MTG6'
            %c = 'kgggggbbbbbrrrrr';
            c = [f f f;
                 f/2 1 f/2;
                 f 1 f;
                 f 1 f;
                 f 1 f;
                 f 1 f;
                 f/2 f/2 1;
                 f f 1;
                 f f 1;
                 f f 1;
                 f f 1;
                 1 f/2 f/2;
                 1 f f;
                 1 f f;
                 1 f f
                 1 f f];
        case 'FLL3'
            %c = 'kgggggbbbbrrrrr';
            c = [f f f;
                 f/2 1 f/2;
                 f 1 f;
                 f 1 f;
                 f 1 f;
                 f 1 f;
                 f/2 f/2 1;
                 f f 1;
                 f f 1;
                 f f 1;
                 1 f/2 f/2;
                 1 f f;
                 1 f f;
                 1 f f
                 1 f f];
        case 'MFG3'
            %c = 'kgggggbbbbbrrrrr';
            c = [f f f;
                 f/2 1 f/2;
                 f 1 f;
                 f 1 f;
                 f 1 f;
                 f 1 f;
                 f/2 f/2 1;
                 f f 1;
                 f f 1;
                 f f 1;
                 f f 1;
                 1 f/2 f/2;
                 1 f f;
                 1 f f;
                 1 f f
                 1 f f];
        otherwise
            disp('Invalid type of cell for color coding');
            quit
    end
    c = flipud(c);
    hold on
    boxplot(x.*sf,gr,'Symbol','');
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),c(j,:),'FaceAlpha',.5);
    end
    ylim([1e-12 1e-3].*sf);
    yticks([1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3]);
    set(gca, 'YScale', 'log');
    posref = get(gca,'position');
    maxh = 7.488000000000001e+02;
    set(gcf,'position',[posref(1) posref(2) maxh 2.864275440684277e+02]);
    set(gcf,'color','w');
    view([90 90])
    xtickangle(0);
    title(groups{g+1});
    saveas(fig,[groups{g+1} '_CSP.svg'],'svg');
    savefig([groups{g+1} '_CSP']);
    hold off
end
posref = get(gcf,'position');
posgca = get(gca,'position');
close all
gaat1 = cell(1,length(allagt1));
gaat1(:) = {['Aggregate (' types{1} ')']};
gsmt1 = cell(1,length(somagt1));
gsmt1(:) = {['Soma (' types{1} ')']};
gamt1 = cell(1,length(axonagt1));
gamt1(:) = {['Axon (' types{1} ')']};
gdmt1 = cell(1,length(dendagt1));
gdmt1(:) = {['Dendrites (' types{1} ')']};

gaat2 = cell(1,length(allagt2));
gaat2(:) = {['Aggregate (' types{2} ')']};
gsmt2 = cell(1,length(somagt2));
gsmt2(:) = {['Soma (' types{2} ')']};
gamt2 = cell(1,length(axonagt2));
gamt2(:) = {['Axon (' types{2} ')']};
gdmt2 = cell(1,length(dendagt2));
gdmt2(:) = {['Dendrites (' types{2} ')']};

gaat3 = cell(1,length(allagt3));
gaat3(:) = {['Aggregate (' types{3} ')']};
gsmt3 = cell(1,length(somagt3));
gsmt3(:) = {['Soma (' types{3} ')']};
gamt3 = cell(1,length(axonagt1));
gamt3(:) = {['Axon (' types{3} ')']};
gdmt3 = cell(1,length(dendagt3));
gdmt3(:) = {['Dendrites (' types{3} ')']};

gaat4 = cell(1,length(allagt4));
gaat4(:) = {['Aggregate (' types{4} ')']};
gsmt4 = cell(1,length(somagt4));
gsmt4(:) = {['Soma (' types{4} ')']};
gamt4 = cell(1,length(axonagt4));
gamt4(:) = {['Axon (' types{4} ')']};
gdmt4 = cell(1,length(dendagt4));
gdmt4(:) = {['Dendrites (' types{4} ')']};

merge = {gaat1,gaat2,gaat3,gaat4,gsmt1,gsmt2,gsmt3,gsmt4,gamt1,...
    gamt2,gamt3,gamt4,gdmt1,gdmt2,gdmt3,gdmt4};
gr = horzcat(merge{:});

f = 0.7;
t = 1;
if compareABSVAL
    x = cat(2,abs(allagt1),abs(somagt1),abs(axonagt1),abs(dendagt1),...
        abs(allagt2),abs(somagt2),abs(axonagt2),abs(dendagt2),...
        abs(allagt3),abs(somagt3),abs(axonagt3),abs(dendagt3),...
        abs(allagt4),abs(somagt4),abs(axonagt4),abs(dendagt4));
else
    x = cat(2,allagt1,allagt2,allagt3,allagt4,...
        somagt1,somagt2,somagt3,somagt4,...
        axonagt1,axonagt2,axonagt3,axonagt4,...
        dendagt1,dendagt2,dendagt3,dendagt4);
end
fig = figure;
c = [f f f;
     f f f;
     f f f;
     f f f;
     f/2 1 f/2;
     f/2 1 f/2;
     f/2 1 f/2;
     f/2 1 f/2;
     f/2 f/2 1;
     f/2 f/2 1;
     f/2 f/2 1;
     f/2 f/2 1;
     1 f/2 f/2;
     1 f/2 f/2;
     1 f/2 f/2;
     1 f/2 f/2];
c = flipud(c);
boxplot(x.*sf,gr,'Symbol','');
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),c(j,:),'FaceAlpha',.5);
end
set(gca, 'YScale', 'log');
set(gcf,'color','w');

ylim([1e-12 1e-3].*sf);
yticks([1e-12 1e-11 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3].*sf);
set(gcf,'position',[posref(1) posref(2) maxh 2.864275440684277e+02]);
view([90 90])
xtickangle(0);
title('All Analyzed Cell Types');
saveas(fig,'AllTypes_CSP.svg','svg');
savefig('AllTypes_CSP.fig');
hold off

% Perform ANOVA testing and output results
xa = cat(1,[gamt1,gamt2,gamt3,gamt4]);
aa = cat(1,axonagt1',axonagt2',axonagt3',axonagt4');
pa = anova1(abs(aa)*sf,xa);

xg = cat(1,[gaat1,gaat2,gaat3,gaat4]);
ag = cat(1,allagt1',allagt2',allagt3',allagt4');
pg = anova1(abs(ag)*sf,xg);

xd = cat(1,[gdmt1,gdmt2,gdmt3,gdmt4]);
ad = cat(1,dendagt1',dendagt2',dendagt3',dendagt4');
pd = anova1(abs(ad)*sf,xd);

xs = cat(1,[gsmt1,gsmt2,gsmt3,gsmt4]);
as = cat(1,somagt1',somagt2',somagt3',somagt4');
ps = anova1(abs(as)*sf,xs);

disp('soma:')
disp(ps)
disp('axon:')
disp(pa)
disp('agg:')
disp(pg)
disp('dend:')
disp(pd)
