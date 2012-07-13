%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  File:  PlotBops.m (beta)                                 %
%  Author: Anupam Karmakar, a.karmakar@fz-juelich.de        %
%  This script plots xxxx.xy file sequence and creates mpeg %
%  encoded movie bopsmovie_xxxx.avi.                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, clf;clc
figure(1)
%%
% snapshots to plot: start, end
start_timesequence=01; % Start no. of time sequence
end_timesequence=30; % End no. of time sequence
%
name='eysi'; % Which series to plot? e.g. exsi, eysi, nenc, ninc 
% axis set to autoscale or fix axis ranges.
axes='auto';% Manual/auto axes
%axes='manual';% Manual/auto axes
% If manual axis set!
xmin=0; % x-axis min value
xmax=30;% x-axis max value
%
ymin=-10.; % y-axis min value
ymax=+10; % y-axis max value
% Now we have enough info:
%%
if start_timesequence == 0 %some xxxx00.xy cases
    get00 ='yes'; % Does filename xxx000.xy exists?
else get00 = 'no';
end
offset=start_timesequence -1;
numframes=end_timesequence - start_timesequence; % no. of timeframes.
A=moviein(numframes); % create the movie matrix
set(gcf,'NextPlot','replacechildren')
axis auto % fix the axes

%%
%file iterations
%
for i=start_timesequence:end_timesequence
    fileindex=int2str(i);
        if strcmp(get00,'yes')
            myfilename = strcat(name,num2str(i-1,'%02i'),'.xy');
    else
        myfilename = strcat(name,num2str(i,'%02i'),'.xy');
    end
    mydata = load (myfilename);     %  read data into plotBOPS matrix
    x = mydata(:,1);                %  copy first column
    field=mydata(:,2);              %  and second column 
    if strncmp(name,'pxx',3);
       plot(x,field,'.','markerSize',4);        %  plot in dots if pxxe or pxxi
    elseif strncmp(name,'f',1)
       semilogy(x,field,'-','linewidth',3)      % make ylogscale for energy spectrum
    else
       plot(x,field,'-','linewidth',3);         %  plot field with lines
    end
    if strncmp(name,'f',1)			% make ylogscale for energy spectrum
       xlabel('Energy in MeV','fontsize',16);
    else xlabel('X in \mu','fontsize',16); 	%  add axis labels and plot title
    end
    ylabel(name,'fontsize',16);
%
    if strcmp(axes, 'manual')
    xlim([xmin xmax]);
    ylim([ymin ymax]);
    else axis auto;
    end
%
    grid on;
    title(myfilename,'fontsize',16);
    A(:,i-offset)=getframe(gcf);
    hold off;
end
%%
%movie creation part
%%%%%%%%%%%%%%%%%%%%%%%
movie(A,1,1); % Play the MATLAB movie
movie2avi(A,['bopsmovie_',name,'.avi'],'fps',1);
close all;
%save (['bopsmovie_',name,'.mat'],A); % save the MATLAB movie to a file

% Notice the MPEG file is about a quarter of the size of the MATLAB movie file
implay (['bopsmovie_',name,'.avi']); % Play the MPEG movie 
%%
