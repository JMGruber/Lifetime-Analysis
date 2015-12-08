allfiles =20:33;
readdir = 'O:\Michael\2015\PSII supercomplex Pengqi\Grana from Ricard\experiment +GOC and FC\';
writedir = readdir;
APD = 1; % "1" is MPD, "0" is SPCM AQRH 16
filename2 = 'trace';
timeres = 0.01;
x=.000:.004:13.15;

%subtrace=1;


    

        
    


HistX=0:100:13000000;

A=4*1000024.00*250/250;
B=4*1000024.00*250/333;
%%
C=4*1000022.00*250/433;
D=4*1000024.00*250/475;
G=10000.2454*3;
H=1/30007.489*1.000024E9;
I=1/50000.0153*1.000024E9;

Shift=169.6E5;

number=2;


histdecayX=x(1:length(x)-1).';


parfor tr=allfiles
    for subtrace=6
        
        f=[1 2 2 2 2 150 1 10]; %first is 1000
        Ubound=[1E6 0.2E6 4E4 4E6 2.5E6 35000 20E6]; %first is 5E5
        %XBounds=[0 Ubound(subtrace)];
        F=1/f(subtrace)*0.99772301E9; % 30 kHz:30007.4845 30007.4916 or 1/f(subtrace)*1.00002442E9
        %trial 76E6 MHz F=1/f(subtrace)*0.99784158E9
        %1/f(subtrace)*0.99784159E9
        %10 kHz:0.99805336E9
        %10 kHz for annihilation trace 11 150 nW:0.99919358E9
        %AOM 215 SMA:0.99825305E9
        %30 kHz annihilation PSII 150 nW F=1/f(subtrace)*0.99757519E9;
       %F=1/f(subtrace)*0.9978299E9; 
        
    if number > 1
        add=['_' int2str(subtrace)];
    else
        add='';
    end 
    [trace,delaytimes,RealTimeHist] = read_pt3_v4_AOM(timeres,readdir,[filename2 int2str(tr) add  '.pt3'],APD,F,HistX,XBounds,Shift,subtrace);    %int2str(subtrace)
    
    histdecayY=hist(delaytimes(delaytimes>0),x); 
    histdecayY=histdecayY(1:length(histdecayY)-1).';
    histdecay=[histdecayX histdecayY];

        
    %dlmwrite(fullfile('D:\Temp data','Decay_2'),histdecay,' ');
    dlmwrite(fullfile(writedir,['trace' int2str(tr) add]),trace,' ');
    dlmwrite(fullfile(writedir,['decaytrace' int2str(tr) add]),histdecay,' ');
    dlmwrite(fullfile(writedir,['Delaytimes' int2str(tr)  add]),delaytimes',' ');
    end
    tr
end


