%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  File:  PlotBops.m (beta)                                 %
%  Author: Anupam Karmakar, a.karmakar@fz-juelich.de        %
%  This script plots xxxx.xy file sequence and creates mpeg %
%  encoded movie bopsmovie_xx_xx.avi.                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc;
close all;
% figure(1)
%% RUN DATA INPUT BLOCK
wkdir = '/home/anupam/bops/bops_trunk/foil_fs';   % full path of run dir

% GRAPHICS INPUT:get filename end numbers etc.
series1 = 'eysi';   % Which series to plot?
axes1 = 'auto';     % Manual/auto axes
xmin1 = 0.;         % x-axis min value
xmax1 = 30.;        % x-axis max value
ymin1 = -5.;        % y-axis min value
ymax1 = +5.;        % y-axis max value
%
series2 = 'nenc';     % Which series to plot?
axes2 = 'auto';     % Manual/auto axes
xmin2 = 0;          % x-axis min value
xmax2 = 30;         % x-axis max value
ymin2 = 0.;         % y-axis min value
ymax2 = 700;        % y-axis max value
%
series3 = 'fuep';     % Which series to plot?
axes3 = 'auto';     % Manual/auto axes
xmin3 = 0;          % x-axis min value
xmax3 = 1;          % x-axis max value
ymin3 = 1;          % y-axis min value
ymax3 = 1e6;        % y-axis max value
%
series4 = 'pxxe';     % Which series to plot?
axes4 = 'auto';     % Manual/auto axes
xmin4 = 0;          % x-axis min value
xmax4 = 1;          % x-axis max value
ymin4 = 1;          % y-axis min value
ymax4 = 1e6;        % y-axis max value

%% CONTROLS
fprintf('Current working dir %s \n',wkdir);
readname1 = strcat(series1,'*.xy');readname2 = strcat(series2,'*.xy');
readname3 = strcat(series3,'*.xy');readname4 = strcat(series4,'*.xy');
%
Filelist1 = dir(strcat(wkdir,'/',readname1)); 
Filelist2 = dir(strcat(wkdir,'/',readname2));
Filelist3 = dir(strcat(wkdir,'/',readname3));
Filelist4 = dir(strcat(wkdir,'/',readname4));
%
numframes = length(Filelist1);
assert(numframes >= 1, 'Folder not found or contains no readable data!');

%% SET SCREEN
scrsz = get(0,'ScreenSize');
h=figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2],...
    'Name','1D BOPS Simulation Plots','NumberTitle','off');
set(0,'CurrentFigure',h);
set(gcf,'NextPlot','replacechildren')
%
A=moviein(numframes);                   % create the movie matrix
axis auto                               % fix the axes
for i = 1:length(Filelist1)
    %% LOAD DATA
    mydata1 = load (strcat(wkdir,'/',Filelist1(i).name)); %  read data into plotBOPS matrix
    mydata2 = load (strcat(wkdir,'/',Filelist2(i).name)); %  read data into plotBOPS matrix
    mydata3 = load (strcat(wkdir,'/',Filelist3(i).name)); %  read data into plotBOPS matrix
    mydata4 = load (strcat(wkdir,'/',Filelist4(i).name)); %  read data into plotBOPS matrix
%    
    x1 = mydata1(:,1);                %  copy first column
    field1=mydata1(:,2);              %  and second column
    x2 = mydata2(:,1);                %  copy first column
    field2=mydata2(:,2);              %  and second column
    x3 = mydata3(:,1);                %  copy first column
    field3=mydata3(:,2);              %  and second column
    x4 = mydata4(:,1);                %  copy first column
    field4=mydata4(:,2);              %  and second column  
    
    %% PLOTS
    subplot(2,2,1)
    %    get plot and axes types
    if strncmp(series1,'pxx',3);
        plot(x1,field1,'.','markerSize',4);           %  plot in dots if pxxe or pxxi
    elseif strncmp(series1,'f',1)
        semilogy(x1,field1,'-','linewidth',3)
    else
        plot(x1,field1,'-','linewidth',2);           %  plot field with lines
    end
    if strncmp(series1,'f',1)
        xlabel('Energy in MeV','fontsize',16);
    end
    xlabel('X in \mu','fontsize',16); %  add axis labels and plot title
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
    if strncmp(series2,'pxx',3);
        plot(x2,field2,'.','markerSize',4);           %  plot in dots if pxxe or pxxi
    elseif strncmp(series2,'f',1)
        semilogy(x2,field2,'-','linewidth',3)
    else
        plot(x2,field2,'-','linewidth',2);           %  plot field with lines
    end
    if strncmp(series2,'f',1)
        xlabel('Energy in MeV','fontsize',16);
    end
    xlabel('X in \mu','fontsize',16); %  add axis labels and plot title
    ylabel(series2,'fontsize',16);
    xlabel('X in \mu','fontsize',16); %  add axis labels and plot title
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
    if strncmp(series3,'pxx',3);
        plot(x3,field3,'.','markerSize',4);           %  plot in dots if pxxe or pxxi
    elseif strncmp(series3,'f',1)
        semilogy(x3/12,field3,'-','linewidth',3)
    else
        plot(x2,field3,'-','linewidth',2);           %  plot field with lines
    end
    if strncmp(series3,'f',1)
        xlabel('Energy in MeV','fontsize',16);
    end
    xlabel('X in \mu','fontsize',16); %  add axis labels and plot title
    ylabel(series3,'fontsize',16);
    xlabel('X in \mu','fontsize',16); %  add axis labels and plot title
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
    if strncmp(series4,'pxx',3);
        plot(x4,field4,'.','markerSize',4);           %  plot in dots if pxxe or pxxi
    elseif strncmp(series4,'f',1)
        semilogy(x4,field4,'-','linewidth',3)
    else
        plot(x4,field4,'-','linewidth',2);           %  plot field with lines
    end
    if strncmp(series4,'f',1)
        xlabel('Energy in MeV','fontsize',16);
    end
    xlabel('X in \mu','fontsize',16); %  add axis labels and plot title
    ylabel(series3,'fontsize',16);
    xlabel('X in \mu','fontsize',16); %  add axis labels and plot title
    ylabel(series4,'fontsize',16);
    if strcmp(axes4, 'manual')
        xlim([xmin4 xmax4]);
        ylim([ymin4 ymax4]);
    else axis auto;
    end
    title(Filelist4(i).name,'fontsize',16);
    grid on;
    % store movie
    %     A(:,i-offset)=getframe(gcf);
    A(:,i)=getframe(gcf);
    hold off;
end
%% movie creation part
movie(A,1,1); % Play the MATLAB movie
movfilename = strcat(wkdir,'/','bopsmovie_',series1,'_',series2,'_',....
        series3,'_',series4,'.avi');
fprintf('writing movie in %s \n',movfilename)
movie2avi(A,movfilename,'fps',1);
close all;
% save (['bopsmovie_',name,'.mat'],A); % save the MATLAB movie to a file
% Notice the MPEG file is about a quarter of the size of the MATLAB movie file
% implay (movfilename); % Play the MPEG movie
%%