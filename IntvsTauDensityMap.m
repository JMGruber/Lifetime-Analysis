% Draws "first Panel" as density map, i.e., intensity level vs. dwell time

% Just run (as IntvsTauDensityMap)
% Run together with dataDensityLog.m

%% parameter initialisation

x0 = .015; % starting x of input data
xinc = log(1.25); % increment of input x-data

mapxmin = .01; % in sec
mapxmax = 100; % in sec
mapymin = 0;   % in counts / time step
mapymax = 100; % in counts / time step
corrbg = 0; % correct background if too much/little has been subtracted (in counts / time step)

fudgefactor = .1; % smearing factor

%%

x=alldwelltime;
y=allint-corrbg;

ymax=round(max(y));
xdata=x0*(exp(0:xinc:9));

dmap1=dataDensityLog(x,y,ymax,xdata,[mapxmin*intbin mapxmax mapymin mapymax],fudgefactor);

deltay=(ymax-min(y))/size(dmap1,1);
ydata=min(y)+deltay:deltay:ymax;

figure1 = figure('Renderer','painters','InvertHardcopy','off','Color',[1 1 1]);

axes1 = axes('Parent',figure1,'XScale','log','XMinorTick','on',...
    'YTick',mapymin:(mapymax-mapymin)/5:mapymax,...   % [0 20 40 60 80 100]
    'TickDir','out',...    %'Position',[0.06 0.65 0.2 0.3],'Position',[0.3 0.65 0.2 0.3],'Position',[0.5 0.65 0.2 0.3],'Position',[0.7 0.65 0.2 0.3],'Position',[0.06 0.3 0.2 0.3],'Position',[0.3 0.29 0.2 0.3],'Position',[0.5 0.29 0.2 0.3],'Position',[0.7 0.29 0.2 0.3]
    'LineWidth',1,'Layer','top',...
    'FontSize',16,...
    'FontName','arial');
xlim(axes1,[mapxmin*intbin mapxmax]);
ylim(axes1,[mapymin ymax]);
box(axes1,'on');
hold(axes1,'all');

contour(xdata,ydata,dmap1,'LineStyle','none','LevelStep',round(max(max(dmap1))/100),'Fill','on','Parent',axes1);

xlabel('Dwell time (s)','FontSize',16,'FontName','arial');
ylabel(['Intensity level (c/',int2str(timeres*1000),' ms)'],'FontSize',16,'FontName','arial');
saveas(figure1,fullfile(writedir,'IntvsTauDensityMap1.jpg'));
saveas(figure1,fullfile(writedir,'IntvsTauDensityMap1.fig'));


%2nd time with larger fudge (smearing) factor
% dmap1=dataDensityLog(x,y,ymax,xdata,[0.01*intbin 100 0 100],0.3);
% deltay=(ymax-min(y))/size(dmap1,1);
% ydata=min(y)+deltay:deltay:ymax;
% h2 = figure('Renderer','painters','InvertHardcopy','off','Color',[1 1 1]);
% axes1 = axes('Parent',h2,'XScale','log','XMinorTick','on',...
%     'YTick',[0 20 40 60 80 100],...
%     'LineWidth',1,'Layer','top',...
%     'FontSize',20,...
%     'FontName','times');
% xlim(axes1,[0.01*intbin 100]);
% ylim(axes1,[0 100]);
% box(axes1,'on');
% hold(axes1,'all');
% contour(xdata,ydata,dmap1,'LineStyle','none','LevelStep',5,'Fill','on');
% xlabel('Dwell time (s)','FontSize',20,'FontName','times');
% ylabel(['Intensity (c/',int2str(timeres*1000),' ms)'],'FontSize',20,'FontName','times');
% saveas(h2,fullfile(writedir,'IntvsTauDensityMap2.jpg'));
% saveas(h2,fullfile(writedir,'IntvsTauDensityMap2.fig'));
