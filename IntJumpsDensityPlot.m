% Draws a map of the intensity before a jump vs. the intensity after a jump for all intensity transitions 
% 1) first a scatter plot
% 2) Then a density map as a pixelled image
% 3) Then a density map with filled contours

% Just run (as IntJumpsDensityPlot)
% Run together with dataDensity.m

%% parameter initialisation

mapxmin = 0;   % in counts / time step
mapxmax = 100; % in counts / time step
mapymin = 0;   % in counts / time step
mapymax = 100; % in counts / time step

fudgefactor = .2; % smearing factor

%%

x=[];
y=[];
ri=1;
for i=1:length(allnc)
    Ii=allint(ri:ri+allnc(i)-2);
    If=allint(ri+1:ri+allnc(i)-1);
%     plot(Ii,If,'.');
%     hold on;
    
    x=[x;Ii];
    y=[y;If];
    
    ri=ri+allnc(i);
end
h=figure;
plot(x,y,'.')
xlabel(['Initial intensity (c/',int2str(timeres*1000*intbin),' ms)'],'FontSize',16,'FontName','times');
ylabel(['Final intensity (c/',int2str(timeres*1000*intbin),' ms)'],'FontSize',16,'FontName','times');
saveas(h,fullfile(writedir,'IntJumps_Scatter.jpg'));
saveas(h,fullfile(writedir,'IntJumps_Scatter.fig'));

maxx=min(max(x),mapxmax);
maxy=min(max(y),mapymax);
dmap=dataDensity(x,y,round(max(x)),round(max(y)),[mapxmin maxx mapymin maxy],fudgefactor);

deltax=(max(x)-min(x))/size(dmap,1);
xdata=min(x)+deltax:deltax:max(x);
deltay=(max(y)-min(y))/size(dmap,1);
ydata=min(y)+deltay:deltay:max(y);

h2 = figure('Renderer','painters','InvertHardcopy','off',...
    'Color',[1 1 1]);
%colormap(jet(round(max(max(dmap)))));


axes1 = axes('Parent',h2,'YTick',mapymin:(maxy-mapymin)/5:maxy,...
    'XTick',mapxmin:(maxx-mapxmin)/5:maxx,...
    'LineWidth',1,...
    'Layer','top',...
    'FontSize',16,...
    'FontName','times');
box(axes1,'on');
hold(axes1,'all');
xlim(axes1,[mapxmin maxx]);
ylim(axes1,[mapymin maxy]);

image(dmap,'Parent',axes1,'CDataMapping','scaled');
xlabel(['Initial intensity (c/',int2str(timeres*1000*intbin),' ms)'],'FontSize',16,'FontName','times');
ylabel(['Final intensity (c/',int2str(timeres*1000*intbin),' ms)'],'FontSize',16,'FontName','times');
saveas(h2,fullfile(writedir,'IntJumps_DensityMap.jpg'));
saveas(h2,fullfile(writedir,'IntJumps_DensityMap.fig'));




h3 = figure('Renderer','painters','InvertHardcopy','off','Color',[1 1 1]);

% Create axes
axes1 = axes('Parent',h3,'LineWidth',1,'Layer','top','FontSize',16,...
    'YTick',mapymin:(maxy-mapymin)/5:maxy,...
    'XTick',mapxmin:(maxx-mapxmin)/5:maxx,...
    'FontName','times');
xlim(axes1,[mapxmin maxx]);
ylim(axes1,[mapymin maxy]);
box(axes1,'on');
hold(axes1,'all');

contour(dmap,'LineStyle','none','LevelStep',2,'Fill','on');
xlabel(['Initial intensity (c/',int2str(timeres*1000*intbin),' ms)'],'FontSize',16,'FontName','times');
ylabel(['Final intensity (c/',int2str(timeres*1000*intbin),' ms)'],'FontSize',16,'FontName','times');
saveas(h3,fullfile(writedir,'IntJumps_DensityMap2.jpg'));
saveas(h3,fullfile(writedir,'IntJumps_DensityMap2.fig'));