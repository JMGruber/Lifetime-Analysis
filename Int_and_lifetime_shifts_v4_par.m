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
% - IntvsTauDensityMap.m -> Draws first Panel as density map
% - IntJumpsDensityPlot.m -> Draws a density map of Int before a jump vs. Int after a jump
% - Calc_TransRates.m; -> Calculates kinetic information, based on a 3-state intensity model
% - IntHist.m -> draws intensity histogram of intensity data points (counts / time resolution)
% - I'll send the file powerlaw.m upon request if you want to do a power-law analysis
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
% - IntvsTauDensityMap.m -> Draws "first Panel" as density map. Change necessary parameters in m file
% - IntJumpsDensityPlot.m -> Draws a density map of Initial Int vs. Final Int. Change necessary parameters in m file.
% - Calc_TransRates.m; -> Calculates kinetic information, based on a 3-state intensity model. Change necessary parameters in m file
% - IntHist2.m -> draws intensity histogram of intensity data points (counts / time resolution). Change necessary parameters in m file.
%
% Tjaart Krüger
% Vrije Universiteit Amsterdam
% (c) 2010-2012
%
% Improvements:
% - option to exclude files
% - add option to average instead of bin and average 1-2-3,2-3-4,3-4-5,...
% - include dark fractions (often identified as photonbursts). photonbursts should be identified on the basis of dwell times + ratio of on:off
% - Don't specify all k's separately but calculate all for 1
% user-defined sigma deviation
% - fprintf: 1...n possible to display?
% - Don't need all temp parameters with variable sizes: calculations cost time.
% - Check testSM and photonburstfactor. First doesn't work well when dim
% fraction is included
% Better to define large vectors/matrices and remove all the zeros in the end!


%% Initialisation
tic;
% profile on
clear all; close all;
allfiles = 424:500; % A:B or matrix [] of individual spec files. Use positive integers only!
skipfiles = []; % ditto
pathname = 'C:\CD\Data WorkingOn\pH 4\new APD\';
readdir=pathname;
%[filename1, pathname]=uigetfile('*.pt3', 'T3 Mode data:', 0, 0);
writedir = 'C:\CD\Results WorkinOn\pH 4\new APD';
%outfile = [pathname filename(1:end-4) '.out'];
filename2 = 'trace';

routerChannel = 1; % 1 for new APD, 2 for old APD

drawlevels = true;  % draw intensity levels
hidelevel = true;	% Hides figure generation for quicker run
drawdensitymaps = true; % draw intensity maps
parallel = true;	% Then uses PARFOR

forceRun = true; % To force code to try and skip all problematic traces and complete run. But runs twice.

ki = [3.25 1.96 1.37 1]; % sigma-deviation for 1..4 consecutive trace points, resp.  Paper: ki=[3.25 1.96 1.37 1] But 10-ms data noisier than "normal" shot noise, prob due to triplets as 1 source
Nf = 2.1; % noise multiplication factor (1.5 - 2 for good S/N)

thrliveI = 700;  % threshold for final intensity level (in cps, background included, intbin excluded)
thrlivet = 10;  % minimal survival time (in seconds)
thr2c = 10000;  % intensity threshold when definitely >1 complex (in cps, background included, intbin excluded)
intbin = 2;   % number of intensity values to be binned together (input)---- for new version the definition is slightly different
timeres = 0.01; % time resolution (excl intbin)
fixbg = true; % usually more accurate
guessbg = 0; % estimated background; used to test for SM (in cps). If fixbg=false, set rather too large than to small
testSM = true; % test for single quantum unit by <=2-time step into quenched state
testphotonburst = true;
photonburstfactor = 10;   % factor by which dwelling time in Q state exceeds dwelling time in unQ states to define photonburst (typically 10)
intermedfactor = 1; % parameter used to test if there are too many unnatural fluctuations (signifying unstable complex)
% decrease if more fluctuations should be incl. (e.g. 0.9)

maxsize = 60/0.01;  % max size of any of the data sets

histbin = 1;  % number of intensity levels to be binned for plots (output)
Qextent = 0.7;  % extent of quenching (for panel 4)

figtitle = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%               Access beyond this point at own risk                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Preallocate matrix sizes
nfiles=length(allfiles);
sortedtraces=zeros(nfiles,5);  % [1..5] = [allusedtraces alldead alldble allphotonburst allusedtracebg]
allnc=zeros(nfiles,1); % number of intensity levels for each complex
alltauon=zeros(maxsize,nfiles);
alltauoff=zeros(maxsize,nfiles);
allinton=zeros(maxsize,nfiles);
%allintoff=ones(maxsize,maxsize)*10*thr2c; % since int can be 0 after bg subtraction
allintoff=ones(maxsize,nfiles)*10*thr2c; % since int can be 0 after bg subtraction
allint=ones(maxsize,nfiles)*10*thr2c;
alldwelltime=zeros(maxsize,nfiles);
allstarttime=zeros(maxsize,nfiles);
AllDelaytimeShiftIndices=zeros(maxsize,nfiles);
DecayX=0.000:0.004:15;
AllIntLevelHist=[];
IntLevelHist=[];
DimHist=zeros(500,length(DecayX),length(allfiles));
Histfiltered2=[];


alltlive=0;

%% Core
i = 1;
skipfiles=sort(skipfiles);
readfile = fullfile(pathname);
timetemp = toc;
if parallel
	
	forceITemp = zeros(length(allfiles),1);
	test3 = cell(length(allfiles),1);
	test4 = zeros(length(allfiles),5);
	test5 = zeros(length(allfiles),1);
	test7 = zeros(length(allfiles),1);
	a = cell(length(allfiles),1);
	test6 = cell(length(allfiles),1);
	ffile = allfiles(1) - 1;
	if allfiles(1) == 0
		ffile = ffile + 1;
	end
	skipfiles = sort(skipfiles);
	parfor tr=allfiles
		try
			complex=tr;
			j=1;
			DelaytimeShiftIndices=[];
			IntLevelHist=[];
			goodfile=true; %i.e. not to be skipped
			% 	 skipfiles=sort(skipfiles); PARFOR
			while (j<=length(skipfiles))&(skipfiles(j)<=complex)&(goodfile)
				if complex==skipfiles(j)
					goodfile=false;
				end
				j=j+1;
			end
			if goodfile
				tempindex = zeros(1,5);
				% 		readfile = fullfile(pathname); PARFOR
				%    trace=dlmread(readfile);
				%    trace=dlmread(strcat(readdir,'\trace (',int2str(i),')'));
				[trace,delaytimes] = read_pt3_v4_channel(timeres,readfile,[filename2 int2str(tr) '.pt3'], routerChannel);
				
				dlmwrite(fullfile(writedir,['trace' int2str(tr)]),trace,' ');
				
				if intbin>1     % bin intensity values
					trace=IntShifts_binints(trace,intbin); % Alternatively, employ IntShifts_avgints.m
				end
				temptimeres=timeres;    %*****************not timeres but switch to real #counts
				trace(:,2)=trace(:,2)*temptimeres;
				guessbgc=guessbg*temptimeres*intbin;
				thrliveIc=thrliveI*temptimeres*intbin;
				
				tlive=IntShifts_TrimTrace(trace,thrliveIc); % trim trace into surviving time
				test7(tr-ffile) = tlive;
				if tlive<thrlivet/trace(1,1)
					tempindex(2)=tr;
				else
					trace=IntShifts_CR(trace,thr2c*temptimeres*intbin,tlive); % remove excessively large intensities
					dble=false; photonburst=false;
					if testSM
						SM=IntShifts_testSM(trace,guessbgc*3,tlive,thrliveIc);
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
						
						% binned delaytimes is calculated before bg is subtracted
						
						for j = 1:length(intstart);
							DelaytimeShiftIndices(1)=1;
							%DelaytimeShiftIndices(j) = sum(trace(intstart(j):intstart(j)+inttimes(j)-1,2));  %these are the sizes, not really the indices - see next command
							DelaytimeShiftIndices(j+1) = sum(trace(intstart(1):intstart(j)+inttimes(j)-1,2));
							if DelaytimeShiftIndices(j)==0
								IntLevelHist = [IntLevelHist; zeros(1,length(DecayX))];
							else
								IntLevelHist = [IntLevelHist; hist(delaytimes(DelaytimeShiftIndices(j):DelaytimeShiftIndices(j+1)),DecayX)];
							end
						end
						%DelaytimeShiftIndices = [1 DelaytimeShiftIndices(1:end-1)+1];
						
						
						
						
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
						% intlevels(intlevels<-2)=0; % can remove this command, if needed
						inttimes = inttimes.*trace(1,1);  %is this simply times timeres??
						intstart = intstart.*trace(1,1);  %is this simply times timeres??
						alltlive = alltlive+tlive;
						
						% calculations for 4th panel of figure
						if ~dble
							% tempalltauon = zeros(maxsize,1);
							% tempalltauoff = zeros(maxsize,1);
							% tempallinton = zeros(maxsize,1);
							% tempallintoff = zeros(maxsize,1);
							% tempallint = ones(maxsize,1)*10*thr2c;
							% tempalldwelltime = zeros(maxsize,1);
							% tempallstarttime = zeros(maxsize,1);
							
							P4vars = IntShifts_P4vars(intlevels,inttimes,thrliveIc-bg,Qextent);
							t2Level = P4vars(:,1);
							t2Times = P4vars(:,2);
							t2Thresh = P4vars(:,3);
							t2Q = P4vars(:,4);
							le = 0;
							if ~isempty(P4vars)
								le = size(P4vars,1);
								test5(tr-ffile) = le;
								%                         tempalltauon(1:l,1) = P4vars(1:end,1);
								%                         alltauon(:,i) = tempalltauon;
								%                         tempalltauoff(1:l,1) = P4vars(1:end,2);
								%                         alltauoff(:,i) = tempalltauoff;
								%                         tempallinton(1:l,1) = P4vars(1:end,3);
								%                         allinton(:,i) = tempallinton;
								%                         tempallintoff(1:l,1) = P4vars(1:end,4);
								%                         allintoff(:,i) = tempallintoff;
								
								% 				  alltauon(1:le,i) = P4vars(:,1);
								% 				  alltauoff(1:le,i) = P4vars(:,2);
								% 				  allinton(1:le,i) = P4vars(:,3);
								% 				  allintoff(1:le,i) = P4vars(:,4);
							end
							
							l = length(intlevels);
							a{tr-ffile} = P4vars(:,1);
							test3{tr-ffile} = {l,P4vars(:,1),P4vars(:,2),P4vars(:,3),...
								P4vars(:,4),intlevels,inttimes,intstart,DelaytimeShiftIndices,...
								IntLevelHist,length(intlevels) tr-ffile};
							%                     tempallint(1:l,1) = intlevels;
							%                     allint(:,i) = tempallint;
							%                     tempalldwelltime(1:l,1) = inttimes;
							%                     alldwelltime(:,i) = tempalldwelltime;
							%                     tempallstarttime(1:l,1) = intstart(1:l);
							%                     allstarttime(:,i) = tempallstarttime;
							% 				allint(1:l,i) = intlevels;
							% 				alldwelltime(1:l,i) = inttimes;
							% 				allstarttime(1:l,i) = intstart;
							% 				AllDelaytimeShiftIndices(1:l+1,i) = DelaytimeShiftIndices;
							% 				AllIntLevelHist=[AllIntLevelHist; IntLevelHist];     PARFOR
							test6{tr-ffile} = {IntLevelHist};
							% 				DimHist(1:size(IntLevelHist,1),1:size(IntLevelHist,2),i)=IntLevelHist;
							dlmwrite(fullfile(writedir,['delaytimes_trace' int2str(tr)]),delaytimes',' ');
							
							
							tempindex(1) = tr;
							tempindex(5) = bg;
							% allnc(i) = length(intlevels);
							
							% allint(1:l,i) = intlevels;       PARFOR
							% alldwelltime(1:l,i) = inttimes;
							% allstarttime(1:l,i) = intstart;
							% AllDelaytimeShiftIndices(1:l,i) = DelaytimeShiftIndices;
							% AllIntLevelHist=[AllIntLevelHist; IntLevelHist];
							% DimHist(1:size(IntLevelHist,1),1:size(IntLevelHist,2),i)=IntLevelHist;
							% dlmwrite(fullfile(writedir,['delaytimes_trace' int2str(tr)]),delaytimes',' ');
							
							% tempindex(1) = tr;
							% tempindex(5) = bg;
							% allnc(i) = length(intlevels);
						end
					end
				end
				test4(tr-ffile,:) = tempindex;
				% 		sortedtraces(i,:) = tempindex;
			end
			fprintf(1,'%s%6.2f\n','core: ',tr);
			% 	 i=i+1;
		catch err
			forceITemp(tr) = tr;
			fprintf('%s\n',['The problem occured during loop number ' int2str(tr)]);
		end
	end
	
	if forceRun	% To force code to try and skip all problematic traces and complete run. But runs twice.
		forceITemp = sort(forceITemp);
		forceI = forceITemp(forceITemp>0);
		if ~isempty(forceI)
			secondRun = true;
			fprintf('\r\n\r\n%s\r\n','Some errors occured. Trying an adapted second run.');
			skipfiles = sort([skipfiles forceI']);
			parfor tr=allfiles
				try
					complex=tr;
					j=1;
					DelaytimeShiftIndices=[];
					IntLevelHist=[];
					goodfile=true; %i.e. not to be skipped
					% 	 skipfiles=sort(skipfiles); PARFOR
					while (j<=length(skipfiles))&(skipfiles(j)<=complex)&(goodfile)
						if complex==skipfiles(j)
							goodfile=false;
						end
						j=j+1;
					end
					if goodfile
						tempindex = zeros(1,5);
						% 		readfile = fullfile(pathname); PARFOR
						%    trace=dlmread(readfile);
						%    trace=dlmread(strcat(readdir,'\trace (',int2str(i),')'));
						[trace,delaytimes] = read_pt3_v4_channel(timeres,readfile,[filename2 int2str(tr) '.pt3'],routerChannel);
						
						dlmwrite(fullfile(writedir,['trace' int2str(tr)]),trace,' ');
						
						if intbin>1     % bin intensity values
							trace=IntShifts_binints(trace,intbin); % Alternatively, employ IntShifts_avgints.m
						end
						temptimeres=timeres;    %*****************not timeres but switch to real #counts
						trace(:,2)=trace(:,2)*temptimeres;
						guessbgc=guessbg*temptimeres*intbin;
						thrliveIc=thrliveI*temptimeres*intbin;
						
						tlive=IntShifts_TrimTrace(trace,thrliveIc); % trim trace into surviving time
						test7(tr-ffile) = tlive;
						if tlive<thrlivet/trace(1,1)
							tempindex(2)=tr;
						else
							trace=IntShifts_CR(trace,thr2c*temptimeres*intbin,tlive); % remove excessively large intensities
							dble=false; photonburst=false;
							if testSM
								SM=IntShifts_testSM(trace,guessbgc*3,tlive,thrliveIc);
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
								
								% binned delaytimes is calculated before bg is subtracted
								
								for j = 1:length(intstart);
									DelaytimeShiftIndices(1)=1;
									%DelaytimeShiftIndices(j) = sum(trace(intstart(j):intstart(j)+inttimes(j)-1,2));  %these are the sizes, not really the indices - see next command
									DelaytimeShiftIndices(j+1) = sum(trace(intstart(1):intstart(j)+inttimes(j)-1,2));
									if DelaytimeShiftIndices(j)==0
										IntLevelHist = [IntLevelHist; zeros(1,length(DecayX))];
									else
										IntLevelHist = [IntLevelHist; hist(delaytimes(DelaytimeShiftIndices(j):DelaytimeShiftIndices(j+1)),DecayX)];
									end
								end
								%DelaytimeShiftIndices = [1 DelaytimeShiftIndices(1:end-1)+1];
								
								
								
								
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
								% intlevels(intlevels<-2)=0; % can remove this command, if needed
								inttimes = inttimes.*trace(1,1);  %is this simply times timeres??
								intstart = intstart.*trace(1,1);  %is this simply times timeres??
								alltlive = alltlive+tlive;
								
								% calculations for 4th panel of figure
								if ~dble
									% tempalltauon = zeros(maxsize,1);
									% tempalltauoff = zeros(maxsize,1);
									% tempallinton = zeros(maxsize,1);
									% tempallintoff = zeros(maxsize,1);
									% tempallint = ones(maxsize,1)*10*thr2c;
									% tempalldwelltime = zeros(maxsize,1);
									% tempallstarttime = zeros(maxsize,1);
									
									P4vars = IntShifts_P4vars(intlevels,inttimes,thrliveIc-bg,Qextent);
									t2Level = P4vars(:,1);
									t2Times = P4vars(:,2);
									t2Thresh = P4vars(:,3);
									t2Q = P4vars(:,4);
									le = 0;
									if ~isempty(P4vars)
										le = size(P4vars,1);
										test5(tr-ffile) = le;
										%                         tempalltauon(1:l,1) = P4vars(1:end,1);
										%                         alltauon(:,i) = tempalltauon;
										%                         tempalltauoff(1:l,1) = P4vars(1:end,2);
										%                         alltauoff(:,i) = tempalltauoff;
										%                         tempallinton(1:l,1) = P4vars(1:end,3);
										%                         allinton(:,i) = tempallinton;
										%                         tempallintoff(1:l,1) = P4vars(1:end,4);
										%                         allintoff(:,i) = tempallintoff;
										
										% 				  alltauon(1:le,i) = P4vars(:,1);
										% 				  alltauoff(1:le,i) = P4vars(:,2);
										% 				  allinton(1:le,i) = P4vars(:,3);
										% 				  allintoff(1:le,i) = P4vars(:,4);
									end
									
									l = length(intlevels);
									a{tr-ffile} = P4vars(:,1);
									test3{tr-ffile} = {l,P4vars(:,1),P4vars(:,2),P4vars(:,3),...
										P4vars(:,4),intlevels,inttimes,intstart,DelaytimeShiftIndices,...
										IntLevelHist,length(intlevels) tr-ffile};
									%                     tempallint(1:l,1) = intlevels;
									%                     allint(:,i) = tempallint;
									%                     tempalldwelltime(1:l,1) = inttimes;
									%                     alldwelltime(:,i) = tempalldwelltime;
									%                     tempallstarttime(1:l,1) = intstart(1:l);
									%                     allstarttime(:,i) = tempallstarttime;
									% 				allint(1:l,i) = intlevels;
									% 				alldwelltime(1:l,i) = inttimes;
									% 				allstarttime(1:l,i) = intstart;
									% 				AllDelaytimeShiftIndices(1:l+1,i) = DelaytimeShiftIndices;
									% 				AllIntLevelHist=[AllIntLevelHist; IntLevelHist];     PARFOR
									test6{tr-ffile} = {IntLevelHist};
									% 				DimHist(1:size(IntLevelHist,1),1:size(IntLevelHist,2),i)=IntLevelHist;
									dlmwrite(fullfile(writedir,['delaytimes_trace' int2str(tr)]),delaytimes',' ');
									
									
									tempindex(1) = tr;
									tempindex(5) = bg;
									% allnc(i) = length(intlevels);
									
									% allint(1:l,i) = intlevels;       PARFOR
									% alldwelltime(1:l,i) = inttimes;
									% allstarttime(1:l,i) = intstart;
									% AllDelaytimeShiftIndices(1:l,i) = DelaytimeShiftIndices;
									% AllIntLevelHist=[AllIntLevelHist; IntLevelHist];
									% DimHist(1:size(IntLevelHist,1),1:size(IntLevelHist,2),i)=IntLevelHist;
									% dlmwrite(fullfile(writedir,['delaytimes_trace' int2str(tr)]),delaytimes',' ');
									
									% tempindex(1) = tr;
									% tempindex(5) = bg;
									% allnc(i) = length(intlevels);
								end
							end
						end
						test4(tr-ffile,:) = tempindex;
						% 		sortedtraces(i,:) = tempindex;
					end
					fprintf(1,'%s%6.2f\n','core: ',tr);
					% 	 i=i+1;
				catch err
					forceITemp(tr) = tr;
					fprintf('%s\n',['The problem occured during loop number ' int2str(tr)]);
				end
			end
		else
			secondRun = false;
		end
	end
	
	skipfiles=sort(skipfiles);
	for i=1:length(test3)
		if ~isempty(test3{i})
			f = cell2mat(test3{i}(12));
			if test5(i) ~= 0
				alltauon(1:test5(i),f) = cell2mat(test3{i}(2));
				alltauoff(1:test5(i),f) = cell2mat(test3{i}(3));
				allinton(1:test5(i),f) = cell2mat(test3{i}(4));
				allintoff(1:test5(i),f) = cell2mat(test3{i}(5));
			end
			l = cell2mat(test3{i}(1));
			allint(1:l,f) = cell2mat(test3{i}(6));
			alldwelltime(1:l,f) = cell2mat(test3{i}(7));
			allstarttime(1:l,f) = cell2mat(test3{i}(8));
			AllDelaytimeShiftIndices(1:l+1,f) = cell2mat(test3{i}(9));
			DimHist(1:size(cell2mat(test3{i}(10)),1),1:size(cell2mat(test3{i}(10)),2),f)=cell2mat(test3{i}(10));
			allnc(f) = cell2mat(test3{i}(11));
			sortedtraces(f,:) = test4(i,:);
			AllIntLevelHist=[AllIntLevelHist; cell2mat(test6{i})];
		end
	end
	alltlive = sum(test7);
	
	% PARALLEL OFF
else
	for tr=allfiles
		complex=tr;
		j=1;
		DelaytimeShiftIndices=[];
		IntLevelHist=[];
		goodfile=true; %i.e. not to be skipped
		skipfiles=sort(skipfiles);
		while (j<=length(skipfiles))&&(skipfiles(j)<=complex)&&(goodfile)
			if complex==skipfiles(j)
				goodfile=false;
			end
			j=j+1;
		end
		if goodfile
			tempindex = zeros(1,5);
			readfile = fullfile(pathname);
			%    trace=dlmread(readfile);
			%    trace=dlmread(strcat(readdir,'\trace (',int2str(i),')'));
			[trace,delaytimes] = read_pt3_v4_channel(timeres,readfile,[filename2 int2str(tr) '.pt3'],routerChannel);
			
			dlmwrite(fullfile(writedir,['trace' int2str(tr)]),trace,' ');
			
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
					SM=IntShifts_testSM(trace,guessbgc*3,tlive,thrliveIc);
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
					
					% binned delaytimes is calculated before bg is subtracted
					
					for j = 1:length(intstart);
						DelaytimeShiftIndices(1)=1;
						%DelaytimeShiftIndices(j) = sum(trace(intstart(j):intstart(j)+inttimes(j)-1,2));  %these are the sizes, not really the indices - see next command
						DelaytimeShiftIndices(j+1) = sum(trace(intstart(1):intstart(j)+inttimes(j)-1,2));
						if DelaytimeShiftIndices(j)==0
							IntLevelHist = [IntLevelHist; zeros(1,length(DecayX))];
						else
							IntLevelHist = [IntLevelHist; hist(delaytimes(DelaytimeShiftIndices(j):DelaytimeShiftIndices(j+1)),DecayX)];
						end
					end
					%DelaytimeShiftIndices = [1 DelaytimeShiftIndices(1:end-1)+1];
					
					
					
					
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
					inttimes = inttimes.*trace(1,1);  %is this simply times timeres??
					intstart = intstart.*trace(1,1);  %is this simply times timeres??
					alltlive = alltlive+tlive;
					
					% calculations for 4th panel of figure
					if ~dble
						%                     tempalltauon = zeros(maxsize,1);
						%                     tempalltauoff = zeros(maxsize,1);
						%                     tempallinton = zeros(maxsize,1);
						%                     tempallintoff = zeros(maxsize,1);
						%                     tempallint = ones(maxsize,1)*10*thr2c;
						%                     tempalldwelltime = zeros(maxsize,1);
						%                     tempallstarttime = zeros(maxsize,1);
						
						P4vars = IntShifts_P4vars(intlevels,inttimes,thrliveIc-bg,Qextent);
						if ~isempty(P4vars)
							l = size(P4vars,1);
							%                         tempalltauon(1:l,1) = P4vars(1:end,1);
							%                         alltauon(:,i) = tempalltauon;
							%                         tempalltauoff(1:l,1) = P4vars(1:end,2);
							%                         alltauoff(:,i) = tempalltauoff;
							%                         tempallinton(1:l,1) = P4vars(1:end,3);
							%                         allinton(:,i) = tempallinton;
							%                         tempallintoff(1:l,1) = P4vars(1:end,4);
							%                         allintoff(:,i) = tempallintoff;
							
							alltauon(1:l,i) = P4vars(:,1);
							alltauoff(1:l,i) = P4vars(:,2);
							allinton(1:l,i) = P4vars(:,3);
							allintoff(1:l,i) = P4vars(:,4);
						end
						
						l = length(intlevels);
						%                     tempallint(1:l,1) = intlevels;
						%                     allint(:,i) = tempallint;
						%                     tempalldwelltime(1:l,1) = inttimes;
						%                     alldwelltime(:,i) = tempalldwelltime;
						%                     tempallstarttime(1:l,1) = intstart(1:l);
						%                     allstarttime(:,i) = tempallstarttime;
						allint(1:l,i) = intlevels;
						alldwelltime(1:l,i) = inttimes;
						allstarttime(1:l,i) = intstart;
						AllDelaytimeShiftIndices(1:l+1,i) = DelaytimeShiftIndices;
						AllIntLevelHist=[AllIntLevelHist; IntLevelHist];
						DimHist(1:size(IntLevelHist,1),1:size(IntLevelHist,2),i)=IntLevelHist;
						dlmwrite(fullfile(writedir,['delaytimes_trace' int2str(tr)]),delaytimes',' ');
						
						
						tempindex(1) = tr;
						tempindex(5) = bg;
						allnc(i) = length(intlevels);
					end
				end
			end
			sortedtraces(i,:) = tempindex;
		end
		fprintf(1,'%s%6.2f\n','core: ',tr);
		i=i+1;
	end
end
partime1 = toc - timetemp;

% remove zeros from matrices
allint2=allint;
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
AllIntLevelHist=[allint AllIntLevelHist(:,1:3291)];
Histindex=find(max(AllIntLevelHist(:,2:end).')>15);
%Histfiltered=AllIntLevelHist(max(AllIntLevelHist(:,2:end).')>15,:); %filter out only the histograms with maximum values greater than 15

[Histfilteredindices(:,1),Histfilteredindices(:,2)]=find(max(DimHist,[],2)>15);
for k=1:length(Histfilteredindices(:,1))
	Histfiltered2(k,:)=DimHist(Histfilteredindices(k,1),:,Histfilteredindices(k,2));
	
end
Histfiltered2=[allint(Histindex) Histfiltered2(:,1:3291)];

%Histfiltered2=[allint Histfiltered2];
allwelltimesselection=alldwelltime(max(AllIntLevelHist(:,2:end).')>15,:);
allstarttimeselection=allstarttime(max(AllIntLevelHist(:,2:end).')>15,:);
DecayX=[0 DecayX(1:3291)];
Histfiltered2=[DecayX;Histfiltered2];
[Histfiltered,Isort]=sortrows(Histfiltered2);

%% Plot output
if ~isempty(allint)
	
	if drawlevels
		timetemp = toc;
		% 	 if parallel
		% 		for i=1:length(allnc)
		% 		  rf=ri+allnc(i)-1;
		% 		  testing(i,1) = i;
		% 		  testing(i,2) = ri;
		% 		  testing(i,3) = rf;
		% 		  ri=ri+allnc(i);
		% 		end
		% 		parfor i=1:length(allnc)
		% 		  %            trace=dlmread(strcat(readfile,[' (' int2str(allusedtraces(i)) ')']));
		% 		  trace=dlmread(strcat(readfile,int2str(allusedtraces(i))));
		% 		  if intbin>1     % bin intensity values
		% 			 trace=IntShifts_binints(trace,intbin); % Alternatively, employ IntShifts_avgints.m
		% 		  end
		% 		  if levelfighide
		% 			 figure('Visible', 'off');
		% 			 fprintf(1,'%s%6.2f\n','draw levels: ',i);
		% 		  else
		% 			 figure;
		% 		  end
		% 		  h=plot(trace(:,1),trace(:,2).*timeres-allusedtracebg(i),'g');
		% 		  hold on;
		% 		  for k = testing(i,2):testing(i,3)
		% 			 x = [allstarttime(k) allstarttime(k)+alldwelltime(k)];
		% 			 y = [allint(k) allint(k)];
		% 			 plot(x,y,'LineWidth',4,'Color','k');
		% 		  end
		% 		  saveas(h,fullfile(writedir,['trace' int2str(allusedtraces(i)) '.jpg']));
		% 		  saveas(h,fullfile(writedir,['trace' int2str(allusedtraces(i)) '.pdf']));
		% 		  if mod(i,10)==0
		% 			 close all;
		% 		  end
		% 		end
		% 	 else
		readfile = fullfile(writedir,filename2); % int2str(i)]);
		ri=1;
		timetemp = toc;
		if parallel
			testing = zeros(length(allnc),2);
			for i = 1:length(allnc)
				rf=ri+allnc(i)-1;
				testing(i,1) = ri;
				testing(i,2) = rf;
				ri=ri+allnc(i);
			end
			parfor tr=1:length(allnc)
				% 		  rf=ri+allnc(tr)-1;
				%            trace=dlmread(strcat(readfile,[' (' int2str(allusedtraces(tr)) ')']));
				trace=dlmread(strcat(readfile,int2str(allusedtraces(tr))));
				if intbin>1     % bin intensity values
					trace=IntShifts_binints(trace,intbin); % Alternatively, employ IntShifts_avgints.m
				end
				if hidelevel || parallel
					figure('Visible','off');
					fprintf(1,'%s%6.2f\n','draw levels: ',tr);
				else
					figure;
				end
				h=plot(trace(:,1),trace(:,2)-allusedtracebg(tr),'g');
				hold on;
				for k = testing(tr,1):testing(tr,2)
					x = [allstarttime(k) allstarttime(k)+alldwelltime(k)];
					y = [allint(k)*100 allint(k)*100];
					plot(x,y,'LineWidth',4,'Color','k');
				end
				ri=ri+allnc(tr);
				saveas(h,fullfile(writedir,['trace' int2str(allusedtraces(tr)) '.jpg']));
				saveas(h,fullfile(writedir,['trace' int2str(allusedtraces(tr)) '.pdf']));
				if mod(tr,10)==0
					close all;
				end
			end
		else
			for tr=1:length(allnc)
				rf=ri+allnc(tr)-1;
				%            trace=dlmread(strcat(readfile,[' (' int2str(allusedtraces(tr)) ')']));
				trace=dlmread(strcat(readfile,int2str(allusedtraces(tr))));
				if intbin>1     % bin intensity values
					trace=IntShifts_binints(trace,intbin); % Alternatively, employ IntShifts_avgints.m
				end
				if hidelevel
					figure('Visible','off');
					fprintf(1,'%s%6.2f\n','draw levels: ',tr);
				else
					figure;
				end
				h=plot(trace(:,1),trace(:,2)-allusedtracebg(tr),'g');
				hold on;
				for k = ri:rf
					x = [allstarttime(k) allstarttime(k)+alldwelltime(k)];
					y = [allint(k)*100 allint(k)*100];
					plot(x,y,'LineWidth',4,'Color','k');
				end
				ri=ri+allnc(tr);
				saveas(h,fullfile(writedir,['trace' int2str(allusedtraces(tr)) '.jpg']));
				saveas(h,fullfile(writedir,['trace' int2str(allusedtraces(tr)) '.pdf']));
				if mod(tr,10)==0
					close all;
				end
			end
		end
	end
	partime2 = toc - timetemp;
end

% Global overview of switches
h=figure; subplot(2,2,1); semilogx(alldwelltime,allint,'.');
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
for tr=1:length(allint)
	p=ceil((allint(tr)-minbin)./histbin);
	if (p<=0)||(isnan(p)), p=1; end
	%        if p>maxbin, p=maxbin; end
	bin.times(p)=bin.times(p)+alldwelltime(tr);
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

if drawdensitymaps
	IntvsTauDensityMap; % First Panel as density map
	IntJumpsDensityPlot; % Draws a density map of Initial Int vs. Final Int
end

Calc_TransRates;

%IntHist2;  %need to change readdir!!



% end

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
	
	dlmwrite(strcat(writedir,'\AllDelaytimeShiftIndices'),AllDelaytimeShiftIndices',' ');
	
	dlmwrite(strcat(writedir,'\complexes used'),allusedtraces','');
	dlmwrite(strcat(writedir,'\complexes dead'),alldead','');
	dlmwrite(strcat(writedir,'\complexes double'),alldble','');
	dlmwrite(strcat(writedir,'\Decay histograms for intensity levels.csv'),Histfiltered','\t');
	
	save(strcat(writedir,'\analysis.mat'));
	
	log = fopen(strcat(writedir,'\computation log.txt'),'a+');
	if parallel
		fprintf(log,'%s\r\n','Parallel on');
	else
		fprintf(log,'%s\r\n','Parallel off');
	end
	fprintf(log,'%s\t%6.3f\t\t%s\t%6.3f\r\n%s\t%6.3f\t\t%s%6.3f\r\n','Core run time:',...
		partime1,'Relative core run time:',nfiles/partime1,'Level run time:',...
		partime2,'Relative level run time:',length(allnc)/partime2);
	fclose(log);
end
% profile viewer

if forceRun
	if secondRun
		fprintf('%s\n','Adapted second run obviously successful.');
		fprintf('%s\n','The following traces where identified as problematic:');
		for i = 1:length(forceI)
			fprintf('%s%s\n','Trace ',int2str(forceI(i)));
		end
	end
	if secondRun == false
		fprintf('%s\n','Second run not atempted, not problematic trace could be detected.');
	end
end

toc