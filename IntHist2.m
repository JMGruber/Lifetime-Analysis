% intensity histogram of short binned intensities;
% use after execution of IntShifts.m
% Just run (as IntHist2)

% To do: get Ithr from double fit
% fit: http://www.mathworks.com/products/statistics/demos.html?file=/products/demos/shipping/stats/customdist1demo.html

%% parameter initialisation

intcal=1; %intensity calibration (if you're sure the setup was misaligned)
%intbin=1; % bin input data into equal time bins - done in main script
%histbin=1; % bin output data - done in main script
%numbertoavg=10; %normalization as in Vanden Bout et al. 1998; risky to use, since many traces don't start in low-emitting state

%% 

intdata=[];
norms=[];normint=[];


for i=1:length(allusedtraces)
    j=allusedtraces(i);
%    trace=dlmread(strcat(readdir,'\trace (',int2str(i),')'));
    trace=dlmread(strcat(readdir,'\trace',int2str(j)));
    if intbin>1     % bin intensity values
        trace=IntShifts_binints(trace,intbin); % Alternatively, employ IntShifts_avgints.m
    end
    trace(:,2)=(trace(:,2)*timeres-allusedtracebg(i))*intcal;
    thrlive=(thrliveI*timeres*intbin-allusedtracebg(i))*intcal;
%    trace(:,2)=(trace(:,2)*timeres)*intcal;
%    thrlive=thrliveI*timeres*intbin;
    tlive=IntShifts_TrimTrace(trace,thrlive); % trim trace into surviving time
    if tlive>0
        int=trace(1:tlive,2);
        intdata=[intdata int'];
    end
    fprintf(1,'%6.2f\n',j);
end
xx=ceil(min(intdata)):histbin:floor(max(intdata));
freq=hist(intdata,xx);

if intcal>1 %getting rid of zero ints
    for i=2:length(freq)-1
        if freq(i)==0
            freq(i)=(freq(i-1)+freq(i+1))/2;
        end
    end
end

x=min(intdata):max(intdata);
h=figure; hist(intdata,xx); 
xlabel('Counts (per 20 ms)'); ylabel('# events');
saveas(h,fullfile(writedir,'hist.jpg'));
dlmwrite(strcat(writedir,'\','hist results'),[xx;freq]','delimiter',' ');