% This script resolves the intensity levels of intensity-time traces, 
% performs various calculations and gives different figures as output
% 
% Intensities are resolved by using the simple algorithm as described in:
% Krüger, Ilioaia, and Van Grondelle,
% Fluorescence Intermittency from the Main Plant Light-Harvesting Complex: Resolving Shifts between Intensity Levels
% Journal of Physical Chemistry B 115:5071–5082, 2011.
%
% Additions to previous version:
% - option to draw all int levels of all traces
% - option to fix background - often more accurate
% - IntvsTauContourPlot.m -> Draws first Panel as density map
% - IntJumpsDensityPlot.m -> Draws a density map of Initial Int vs. Final Int
% - Calc_TransRates.m; -> Calculates kinetic information, based on a 3-state intensity model
% - IntHist.m -> draws intensity histogram of intensity data points (counts / time resolution)
% - I'll send the file powerlaw.m on request if you want to do a power-law analysis
%
% Notes: very last int level is not resolved, since a nonblinking event cuts it short (either end of acquisition time or bleaching)
%
% First execute 'matlabpool' to activate parallel computing
% Code written for Parallel Computing: n CPUs will decrease execution time by a factor n.
% (though I have not checked if all permutations work for parallel computing)
%
% Run together with:
% - IntShifts_Algorithm.m
% - IntShifts_binints.m
% - IntShifts_CR.m
% - IntShifts_P4vars.m
% - IntShifts_TestSM.m
% - IntShifts_TrimTrace.m
% - IntvsTauContourPlot.m -> Draws first Panel as density map
% - IntJumpsDensityPlot.m -> Draws a density map of Initial Int vs. Final Int
% - Calc_TransRates.m; -> Calculates kinetic information, based on a 3-state intensity model
% - IntHist.m -> draws intensity histogram of intensity data points (counts / time resolution)
% 
% Tjaart Krüger
% Vrije Universiteit Amsterdam
% (c) 2010-2012
%
% Improvements: 
% - option to exclude files
% - include dark fractions (often identified as photonbursts). photonbursts should be identified on the basis of dwell times + ratio of on:off
% - Don't specify all k's separately but calculate all for 1
% user-defined sigma deviation
% - fprintf: 1...n possible to display?
% - Don't need all temp parameters with variable sizes: calculations cost time.
% Better to define large vectors/matrices and remove all the zeros in the end!


%% Initialisation
tic;
clear all; close all;
startfile=5; endfile=304; % files must be numbered consecutively!
skipfiles=[];
includefiles=[175:186];
readdir='J:\2014\Joshua and Herman\pH dependence PCD PCA\pH 9';
writedir='J:\2014\Joshua and Herman\pH dependence PCD PCA\pH 9\trace analysis';
filename='trace';

ki = [3.25 1.96 1.37 1]; % sigma-deviation for 1..4 consecutive trace points, resp.  Paper: ki=[3.25 1.96 1.37 1] But 10-ms data noisier than "normal" shot noise, prob due to triplets as 1 source
Nf = 2.5; % noise multiplication factor (1.0 - 1.1 for good SNR)

thrliveI = 4000;  % threshold for final intensity level (in cps, background included, intbin excl)
thrlivet = 5;  % minimal survival time (in seconds)
thr2c = 10000;  % intensity threshold when definitely >1 complex (in cps, background included, intbin excl)
intbin = 1;   % number of intensity values to be binned together (input)
timeres = 0.01; % time resolution (excl intbin)
fixbg = true; % usually more accurate
guessbg = 100; % estimated background; used to test for SM (in cps). If fixbg=false, set rather too large than to small
testSM = true; % test for single quantum unit by <=2-time step into quenched state
testphotonburst = false;
photonburstfactor = 10;   % factor by which dwelling time in Q state exceeds dwelling time in unQ states to define photonburst (typically 10)
intermedfactor = 1; % parameter used to test if there are too many unnatural fluctuations (signifying unstable complex)
                    % decrease if more fluctuations should be incl. (e.g. 0.9)

maxsize = 60/0.01;  % max size of any of the data sets

histbin = 1;  % number of intensity levels to be binned for plots (output)
Qextent = 0.7;  % extent of quenching (for panel 4)

figtitle = '';

drawlevels = 1;  % draw intensity levels of the first indexed complex

% changed trace=dlmread(readfile,'',2,0); (delete last parts) line120
% trace(:,2)=trace(:,2)*1000; line 126
% same 250 + 251
% change in inhist2.m:
% trace=dlmread(strcat(readdir,'\trace',int2str(j)),'',2,0);
%                            trace(:,2)=trace(:,2)*1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Access beyond this point only if you know what you're doing      %%%                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Preallocate matrix sizes
nfiles=endfile-startfile+1;
sortedtraces=zeros(nfiles,5);  % [1..5] = [allusedtraces alldead alldble allphotonburst allusedtracebg]
allnc=zeros(nfiles,1); % number of intensity levels for each complex
alltauon=zeros(maxsize,nfiles);
alltauoff=zeros(maxsize,nfiles);
allinton=zeros(maxsize,nfiles);
allintoff=ones(maxsize,maxsize)*10*thr2c; % since int can be 0 after bg subtraction
allint=ones(maxsize,nfiles)*10*thr2c;
alldwelltime=zeros(maxsize,nfiles);
allstarttime=zeros(maxsize,nfiles);

alltlive=0;


%% Core
parfor i=1:nfiles
    j=1;
    k=1;
        goodfile=false;
        
      while (j<=length(includefiles))&&(includefiles(j)<=startfile+i-1)&&~(goodfile)
      if startfile+i-1==includefiles(j)
          goodfile=true;
      end
      j=j+1;
      end
    
    while (k<=length(skipfiles))&&(skipfiles(k)<=startfile+i-1)&&(goodfile)
        if startfile+i-1==skipfiles(k)
            goodfile=false;
        end
        k=k+1;
    end
    

    
    if goodfile
    tr = startfile+i-1;
    tempindex = zeros(1,5);
    readfile = fullfile(readdir,[filename int2str(tr)]);
    trace=dlmread(readfile);
%    trace=dlmread(strcat(readdir,'\trace (',int2str(tr),')'));
    if intbin>1     % bin intensity values
        trace=IntShifts_binints(trace,intbin); % Alternatively, employ IntShifts_avgints.m
    end
    temptimeres=timeres;    %*****************not timeres but switch to real #counts
    trace(:,2)=trace(:,2)*temptimeres;
    guessbgc=guessbg*temptimeres*intbin;
    thrliveIc=thrliveI*temptimeres*intbin;
    
    tlive=IntShifts_TrimTrace(trace,thrliveIc); % trim trace into surviving time
    if tlive<thrlivet/trace(1,1);
        tempindex(2)=tr; 
    else
       trace=IntShifts_CR(trace,thr2c*temptimeres*intbin,tlive); % remove excessively large intensities
       dble=false; photonburst=false;
       if testSM
           SM=IntShifts_testSM(trace,guessbgc*5,tlive,thrliveIc);
           if ~SM
              dble=true;
              tempindex(3)=tr;
           end
       end
       if (testphotonburst)&&((length(find(trace(1:tlive,2)<thrliveIc))>photonburstfactor*length(find(trace(1:tlive,2)>thrliveIc)))...
                            ||(sum(trace(1:tlive,2)<intermedfactor*thrliveIc & trace(1:tlive,2)>(guessbgc+thrliveIc)/2)... % intermediate intensities
                            > sum(trace(1:tlive,2)>intermedfactor*thrliveIc | trace(1:tlive,2)<(guessbgc+thrliveIc)/2)))       % are more than Q and unQ
            photonburst=true;
            tempindex(4)=tr; 
       end
        
        if (~photonburst)&&(~dble)
            [intlevels,inttimes,intstart,SM] = IntShifts_Algorithm(trace,ki,Nf,tlive,thr2c*temptimeres*intbin);
            if ~SM
                tempindex(3)=tr; 
                dble=true;
            end

            % estimate background
            if fixbg
                bg=guessbgc;
            else
                if tlive<trace(end,1)/trace(1,1);
                    blI=trace(tlive:end,2);
                    avgblI=sum(blI)./length(blI);
                    [bg,bgi] = min([intlevels avgblI]);
                else
                    [bg,bgi] = min(intlevels);
                end
                if (bgi<=length(inttimes))&&(inttimes(bgi)*intbin<=2)   % if short dwell times don't represent bad estimations, exchange "2" with "1"
                    levels2 = intlevels;
                    times2 = inttimes;
                    while (bgi>1)&&(bgi<length(levels2))&&(length(levels2)>3)&&(times2(bgi)<=2)
                        levels2 = [levels2(1:bgi-1) levels2(bgi+1:length(levels2))];
                        times2 = [times2(1:bgi-1) times2(bgi+1:length(times2))];
                        [bg,bgi] = min(levels2);
                    end
                end
                bg = min(bg,guessbgc);
            end
            intlevels = intlevels-bg;
%            intlevels(intlevels<-2)=0; % can remove this command, if needed
            inttimes = inttimes.*trace(1,1);
            intstart = intstart.*trace(1,1);
            alltlive = alltlive+tlive;
            
            % calculations for 4th panel of figure
            if ~dble
                tempalltauon = zeros(maxsize,1);
                tempalltauoff = zeros(maxsize,1);
                tempallinton = zeros(maxsize,1);
                tempallintoff = zeros(maxsize,1);
                tempallint = ones(maxsize,1)*10*thr2c;
                tempalldwelltime = zeros(maxsize,1);
                tempallstarttime = zeros(maxsize,1);
                
                P4vars = IntShifts_P4vars(intlevels,inttimes,thrliveIc-bg,Qextent);
                if ~isempty(P4vars)
                    l = size(P4vars,1);
                    tempalltauon(1:l,1) = P4vars(1:end,1);     
                    alltauon(:,i) = tempalltauon;
                    tempalltauoff(1:l,1) = P4vars(1:end,2);     
                    alltauoff(:,i) = tempalltauoff;
                    tempallinton(1:l,1) = P4vars(1:end,3);     
                    allinton(:,i) = tempallinton;
                    tempallintoff(1:l,1) = P4vars(1:end,4);     
                    allintoff(:,i) = tempallintoff;
                end
                
                l = length(intlevels);
                tempallint(1:l,1) = intlevels;
                allint(:,i) = tempallint;
                tempalldwelltime(1:l,1) = inttimes;
                alldwelltime(:,i) = tempalldwelltime;
                tempallstarttime(1:l,1) = intstart(1:l);
                allstarttime(:,i) = tempallstarttime;
                
                tempindex(1) = tr; 
                tempindex(5) = bg;
                allnc(i) = length(intlevels);
            end
        end
    end
    sortedtraces(i,:) = tempindex;
    fprintf(1,'%6.2f\n',tr);
    end
end
% remove zeros from matrices
allusedtraces=sortedtraces(sortedtraces(:,1)>0,1);
alldead=sortedtraces(sortedtraces(:,2)>0,2);
alldble=sortedtraces(sortedtraces(:,3)>0,3);
allphotonburst=sortedtraces(sortedtraces(:,4)>0,4);
allusedtracebg=sortedtraces(sortedtraces(:,5)>0,5);
allint=allint(allint<(10*thr2c));
alldwelltime=alldwelltime(alldwelltime>0);
allstarttime=allstarttime(allstarttime>0);
allnc=allnc(allnc>0);
alltauon=alltauon(alltauon>0);
alltauoff=alltauoff(alltauoff>0);
allinton=allinton(allinton>0);
allintoff=allintoff(allintoff<10*thr2c);

%% Plot output
if ~isempty(allint)
    if drawlevels
        readfile = fullfile(readdir,filename); % int2str(tr)]);
        ri=1;
        for i=1:length(allnc)
            rf=ri+allnc(i)-1;
%            trace=dlmread(strcat(readfile,[' (' int2str(allusedtraces(i)) ')']));
            trace=dlmread(strcat(readfile,int2str(allusedtraces(i))));
            if intbin>1     % bin intensity values
                trace=IntShifts_binints(trace,intbin); % Alternatively, employ IntShifts_avgints.m
            end
            figure; h=plot(trace(:,1),trace(:,2).*timeres-allusedtracebg(i),'g');
            hold on;
            for k = ri:rf
                x = [allstarttime(k) allstarttime(k)+alldwelltime(k)];
                y = [allint(k) allint(k)];
                plot(x,y,'LineWidth',4,'Color','k');
            end
            ri=ri+allnc(i);
            saveas(h,fullfile(writedir,['trace' int2str(allusedtraces(i)) '.jpg']));
            saveas(h,fullfile(writedir,['trace' int2str(allusedtraces(i)) '.pdf']));
            if mod(i,10)==0
                close all;
            end
        end
    end

    % Global overview of switches
    h=figure(3); subplot(2,2,1); semilogx(alldwelltime,allint,'.');
    xlabel('Dwell time (s)'); ylabel(['Intensity (c/',int2str(timeres*1000*intbin),' ms)']);
    title(figtitle);

%    histbin=histbin*timeres; % Take note!
    maxbin=ceil(max(allint)/histbin)*histbin;
    %maxbin=ceil(thr2c*timeres/histbin)*histbin;
    minbin=ceil(min(allint)/histbin)*histbin;

    markers=minbin:histbin:maxbin;
    m=length(markers);
    bin.times=zeros(1,m);
    bin.n2=zeros(1,m);
    for i=1:length(allint) 
        p=ceil((allint(i)-minbin)./histbin);
        if (p<=0)||(isnan(p)), p=1; end
%        if p>maxbin, p=maxbin; end
        bin.times(p)=bin.times(p)+alldwelltime(i);
        bin.n2(p)=bin.n2(p)+1;
    end
    bin.times=bin.times./(alltlive*timeres)*100;
    bin.n2=bin.n2./(alltlive*timeres)*60;
    subplot(2,2,2); bar(markers,bin.n2); 
    xlabel(['Intensity (c/',int2str(timeres*1000*intbin),' ms)']); ylabel('Access frequency (min^{-1})');
    subplot(2,2,3); bar(markers,bin.times);
    xlabel(['Intensity (c/',int2str(timeres*1000*intbin),' ms)']); ylabel('Total dwell time (%)');

    edges=timeres*(exp(0:log(2):9));
    tausortoff=histc(alltauoff,edges)./(alltlive*timeres)*60;
    if ~isempty(tausortoff)
        subplot(2,2,4); semilogx(edges,tausortoff,'.-'); ylabel('Switching frequency (min^{-1})'); xlabel('Dwell time (s)');
    end
    tausorton=histc(alltauon,edges)./(alltlive*timeres)*60;
    if ~isempty(tausorton)
        hold on;
        subplot(2,2,4); semilogx(edges,tausorton,'o-');
    end
    saveas(h,fullfile(writedir,'analysis.jpg'));
    
    %IntvsTauContourPlot; % First Panel as density map
    
    %IntJumpsDensityPlot; % Draws a density map of Initial Int vs. Final Int
    
    Calc_TransRates;
    
    IntHist2;
    
end

%% File output
%dlmwrite(strcat(writedir,'\','parameters'),['bg  ',int2str(bg)],'');    
if ~isempty(allusedtraces)
    dlmwrite(strcat(writedir,'\parameters'),['k1 k2 k3 k4 Nf ',int2str(ki(1)),' ',int2str(ki(2)),' ',int2str(ki(3)),' ',int2str(ki(4)),' ',int2str(Nf)],'');
    dlmwrite(strcat(writedir,'\parameters'),['thrliveI thrlivet  ',int2str(thrliveI),' ',int2str(thrlivet)],'-append','delimiter','');
    dlmwrite(strcat(writedir,'\parameters'),['thr2c maxbg ',int2str(thr2c),' ',int2str(guessbg)],'-append','delimiter','');
    dlmwrite(strcat(writedir,'\parameters'),['intbin histbin Qextent  ',int2str(intbin),' ',int2str(histbin),' ',int2str(Qextent)],'-append','delimiter','');

    dlmwrite(strcat(writedir,'\panel1'),[alldwelltime allint],' ');  
    dlmwrite(strcat(writedir,'\panel23'),[markers;bin.n2;bin.times]',' ');
    dlmwrite(strcat(writedir,'\panel4'),[edges' tausortoff tausorton],' ');

    dlmwrite(strcat(writedir,'\complexes used'),allusedtraces',''); 
    dlmwrite(strcat(writedir,'\complexes dead'),alldead',''); 
    dlmwrite(strcat(writedir,'\complexes double'),alldble',''); 

    save(strcat(writedir,'\analysis.mat'));
end
toc