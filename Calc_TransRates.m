% Calculates kinetic information, based on a 3-state model:
% A = unQ state
% B = Q state
% C = intermediate state

% Calculates 
% - transrates = all transition rate constants ()
%   kAB = rate const A->B
%   kBA = rate const B->A
%   etc.
% - switchespersec
% - times = total time spent in each state
% - avgI = average intensity

% doesn't resolve switches within A,B or C - see powerlaw.m how to resolve that

% Determine thresholds with IntvsTauDensityMap1.jpg or dwell time vs resolved int levels

unQmin=55;  % threshold between intermediate and unQ states
Qmax=20; % threshold between intermediate and Q states


kAB=[]; kBA=[]; kAC=[]; kCA=[]; kBC=[]; kCB=[];
ri=1;
allN_AB=0; allN_BA=0; allN_AC=0; allN_CA=0; allN_BC=0; allN_CB=0; 
alltauBA=0; alltauAB=0; alltauAC=0; alltauCA=0; alltauBC=0; alltauCB=0;

for i=1:length(allnc)
    Ii=allint(ri:ri+allnc(i)-1);
    tau=alldwelltime(ri:ri+allnc(i)-1);
    
    swAB=Ii(1:end-1)>=unQmin & Ii(2:end)<Qmax; %all switches A->B
    N_AB=length(swAB(swAB>0)); % #switches A->B
    tauAB=sum(tau(swAB)); % total time A
    
    swBA=Ii(1:end-1)<Qmax & Ii(2:end)>=unQmin; %all switches B->A
    N_BA=length(swBA(swBA>0)); % #switches B->A
    tauBA=sum(tau(swBA)); % total time B
    
   % swAC=Ii(1:end-1)>unQmin & Ii(2:end)>Qmax; %all switches A->C
    swAC=Ii(1:end-1)>=unQmin & Ii(2:end)>Qmax & Ii(2:end)<=unQmin; %all switches A->C
    N_AC=length(swAC(swAC>0)); % #switches A->C
    tauAC=sum(tau(swAC)); % total time A
    
    %swCA=Ii(1:end-1)>Qmax & Ii(2:end)>unQmin; %all switches C->A
    swCA=Ii(1:end-1)>Qmax & Ii(1:end-1)<=unQmin & Ii(2:end)>=unQmin; %all switches C->A
    N_CA=length(swCA(swCA>0)); % #switches C->A
    tauCA=sum(tau(swCA)); % total time C

    %swBC=Ii(1:end-1)<Qmax & Ii(2:end)<unQmin; %all switches B->C
    swBC=Ii(1:end-1)<Qmax & Ii(2:end)>Qmax & Ii(2:end)<=unQmin;
    N_BC=length(swBC(swBC>0)); % #switches B->C
    tauBC=sum(tau(swBC)); % total time B
    
    %swCB=Ii(1:end-1)<unQmin & Ii(2:end)<Qmax; %all switches C->B
    swCB=Ii(1:end-1)>Qmax & Ii(1:end-1)<=unQmin & Ii(2:end)<Qmax; %all switches C->B
    N_CB=length(swCB(swCB>0)); % #switches C->B
    tauCB=sum(tau(swCB)); % total time C
    
    if tauAB>0
        allN_AB=allN_AB+N_AB;
        alltauAB=alltauAB+tauAB;
        kAB=[kAB N_AB/tauAB];
    end
    if tauBA>0
        allN_BA=allN_BA+N_BA;
        alltauBA=alltauBA+tauBA;
        kBA=[kBA N_BA/tauBA];
    end
    if tauAC>0
        allN_AC=allN_AC+N_AC;
        alltauAC=alltauAC+tauAC;
        kAC=[kAC N_AC/tauAC];
    end
    if tauCA>0
        allN_CA=allN_CA+N_CA;
        alltauCA=alltauCA+tauCA;
        kCA=[kCA N_CA/tauCA];
    end
    if tauBC>0
        allN_BC=allN_BC+N_BC;
        alltauBC=alltauBC+tauBC;
        kBC=[kBC N_BC/tauBC];
    end
    if tauCB>0
        allN_CB=allN_CB+N_CB;
        alltauCB=alltauCB+tauCB;
        kCB=[kCB N_CB/tauCB];
    end
    ri=ri+allnc(i);
end

avgkAB=allN_AB/alltauAB;
avgkBA=allN_BA/alltauBA;
avgkAC=allN_AC/alltauAC;
avgkCA=allN_CA/alltauCA;
avgkBC=allN_BC/alltauBC;
avgkCB=allN_CB/alltauCB;

avgkAB=round(avgkAB*100)/100;
avgkBA=round(avgkBA*100)/100;
avgkAC=round(avgkAC*100)/100;
avgkCA=round(avgkCA*100)/100;
avgkBC=round(avgkBC*100)/100;
avgkCB=round(avgkCB*100)/100;

transrates=[avgkAB avgkBA avgkAC avgkCA avgkBC avgkCB];
dlmwrite(strcat(writedir,'\Transition Rates'),['kAB kBA kAC kCA kBC kCB ',int2str(avgkAB),' ',int2str(avgkBA),' ',int2str(avgkAC),' ',int2str(avgkCA),' ',int2str(avgkBC),' ',int2str(avgkCB)],'');

timesA=alldwelltime(allint>=unQmin);
timesB=alldwelltime(allint<=Qmax);
timesC=alldwelltime(allint>Qmax & allint<unQmin);
tottimeA=sum(timesA);
tottimeB=sum(timesB);
tottimeC=sum(timesC);
swA=length(timesA);
swB=length(timesB);
swC=length(timesC);

tottime=sum(alldwelltime);
times=[tottimeA tottimeB tottimeC]./tottime;
switches=[swA swB swC]./(swA+swB+swC);
switchespersec=[swA swB swC]./tottime;

% 2 ways to calculate avg I (2nd is approximation)
avgI=dot(allint,alldwelltime)/tottime;
avgI2= dot(markers,bin.times)/sum(bin.times);