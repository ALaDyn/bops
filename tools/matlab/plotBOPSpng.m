%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  File:  PlotBops.m (beta)                                 %
%  Author: Anupam Karmakar, a.karmakar@fz-juelich.de        %
%  This script plots xxxx.xy file sequences into simple png %
%  images. More subplots can easily be added                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotBOPSpng
clear all; clc
close all;
%% RUN DATA INPUT BLOCK
wkdir='/home/anupam/bops/bops_trunk/robinson1'; % full path of run dir
cleangraphics = 'yes';                          % clear old graphics?
set(0,'defaultfigurevisible','off');             % show /gide plotting window
linewidth = 1;                                  % plot linewidth
%% GRAPHICS INPUT:get filename end numbers etc.
series1='eysi';     % Which series to plot?
axes1='auto';       %'auto'; %manual% Manual/auto axes
xmin1=0.;           % x-axis min value
xmax1=30.;          % x-axis max value
ymin1=-4.;          % y-axis min value
ymax1=+4.;          % y-axis max value
%
series2='nenc';     % Which series to plot?
axes2='auto';       % Manual/auto axes
xmin2=0;            % x-axis min value
xmax2=30;           % x-axis max value
ymin2=0.;           % y-axis min value
ymax2=200;          % y-axis max value
%
series3='fuip';     % Which series to plot?
axes3='auto';       %manual% Manual/auto axes
xmin3=0;            % x-axis min value
xmax3=5;            % x-axis max value
ymin3=1;            % y-axis min value
ymax3=1e6;          % y-axis max value
%
series4='pxxe';     % Which series to plot?
axes4='auto';       %manual% Manual/auto axes
xmin4=0;            % x-axis min value
xmax4=1;            % x-axis max value
ymin4=1;            % y-axis min value
ymax4=1e6;          % y-axis max value
%
%
%% CONTROL
fprintf('Current working dir %s \n',wkdir);
readname1=strcat(series1,'*.xy');readname2=strcat(series2,'*.xy');
readname3=strcat(series3,'*.xy');readname4=strcat(series4,'*.xy');
%
Filelist1 = dir(strcat(wkdir,'/',readname1)); Filelist2 = dir(strcat(wkdir,'/',readname2));
Filelist3 = dir(strcat(wkdir,'/',readname3)); Filelist4 = dir(strcat(wkdir,'/',readname4));

n = length(Filelist1);assert(n >= 1, 'Folder not found or contains no readable data!');
wb = waitbar(0,'Initializing...');
%
if(strcmpi(cleangraphics,'yes'))
    recycle on; delete (strcat(wkdir,'/','bopsFigure*.png'));
end

%% LOAD DATA AND PLOT FIGURES
for i = 1:length(Filelist1)
    mydata1 = load (strcat(wkdir,'/',Filelist1(i).name)); %  read data into plotBOPS matrix
    mydata2 = load (strcat(wkdir,'/',Filelist2(i).name)); %  read data into plotBOPS matrix
    mydata3 = load (strcat(wkdir,'/',Filelist3(i).name)); %  read data into plotBOPS matrix
    mydata4 = load (strcat(wkdir,'/',Filelist4(i).name)); %  read data into plotBOPS matrix
    
    x1 = mydata1(:,1);                %  copy first column
    field1=mydata1(:,2);              %  and second column
    x2 = mydata2(:,1);                %  copy first column
    field2=mydata2(:,2);              %  and second column
    x3 = mydata3(:,1);                %  copy first column
    field3=mydata3(:,2);              %  and second column
    
    x4 = mydata4(:,1);                %  copy first column
    field4=mydata4(:,2);              %  and second column
    %
    scrsz = get(0,'ScreenSize');
    h=figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2],...
        'Name','1D BOPS Simulation Plots','NumberTitle','off');
    set(0,'CurrentFigure',h);
    set(gcf,'NextPlot','replacechildren')
    %
    subplot(2,2,1)
    %    get plot and axes types
    if strncmp(readname1,'pxx',3);
        plot(x1,field1,'.','markerSize',2);             %  plot in dots if pxxe or pxxi
    elseif strncmp(readname1,'f',1)
        semilogy(x1,field1,'-','linewidth',linewidth)
    else
        plot(x1,field1,'-','linewidth',linewidth);      %  plot field with lines
    end
    if strncmp(readname1,'f',1)
        xlabel('Energy in MeV','fontsize',16);
    else
        xlabel('X in \mu','fontsize',16);               %  add axis labels and plot title
    end
    ylabel(series1,'fontsize',16);
    if strcmp(axes1, 'manual')
        xlim([xmin1 xmax1]);
        ylim([ymin1 ymax1]);
    else axis auto;
    end
    title(Filelist1(i).name,'fontsize',16);
    grid on;
    %
    subplot(2,2,2)
    if strncmp(readname2,'pxx',3);
        plot(x2,field2,'.','markerSize',2);             %  plot in dots if pxxe or pxxi
    elseif strncmp(readname2,'f',1)
        semilogy(x2,field2,'-','linewidth',linewidth)
    else
        plot(x2,field2,'-','linewidth',linewidth);      %  plot field with lines
    end
    if strncmp(readname2,'f',1)
        xlabel('Energy in MeV','fontsize',16);
    else
        xlabel('X in \mu','fontsize',16);               %  add axis labels and plot title
    end
    ylabel(series2,'fontsize',16);
    if strcmp(axes2, 'manual')
        xlim([xmin2 xmax2]);
        ylim([ymin2 ymax2]);
    else axis auto;
    end
    title(Filelist2(i).name,'fontsize',16);
    grid on;
    %
    subplot(2,2,3)
    if strncmp(readname3,'pxx',3);
        plot(x3,field3,'.','markerSize',2);             %  plot in dots if pxxe or pxxi
    elseif strncmp(readname3,'f',1)
        semilogy(x3,field3,'-','linewidth',linewidth)
    else
        plot(x2,field3,'-','linewidth',linewidth);      %  plot field with lines
    end
    if strncmp(readname3,'f',1)
        xlabel('Energy in MeV','fontsize',16);
    else
        xlabel('X in \mu','fontsize',16);               %  add axis labels and plot title
    end
    ylabel(series3,'fontsize',16);
    if strcmp(axes3, 'manual')
        xlim([xmin3 xmax3]);
        ylim([ymin3 ymax3]);
    else axis auto;
    end
    title(Filelist3(i).name,'fontsize',16);
    grid on;
    %
    subplot(2,2,4)
    if strncmp(readname4,'pxx',3);
        plot(x4,field4,'.','markerSize',2);             %  plot in dots if pxxe or pxxi
    elseif strncmp(readname4,'f',1)
        semilogy(x4,field4,'-','linewidth',linewidth)
    else
        plot(x4,field4,'-','linewidth',linewidth);      %  plot field with lines
    end
    if strncmp(readname4,'f',1)
        xlabel('Energy in MeV','fontsize',16);
    else
        xlabel('X in \mu','fontsize',16);               %  add axis labels and plot title
    end
    ylabel(series4,'fontsize',16);
    if strcmp(axes4, 'manual')
        xlim([xmin4 xmax4]);
        ylim([ymin4 ymax4]);
    else axis auto;
    end
    title(Filelist4(i).name,'fontsize',16);
    grid on;
    
    %% SAVE PNG FIGURES
    % set(gcf, 'PaperPositionMode', 'auto');
    gfilename = strcat(wkdir,'/','bopsFigure',....
        int2str(i),'.png');
    print('-dpng',gfilename)
    fprintf('writing %s \n',gfilename)
    prog = i/n;
    waitbar(prog,wb,strcat('Progress : ',num2str(uint8(prog*100)),' %'));
    close(gcf);
end
%close the waitbar
close(wb);
close all;clear all;
%%
end