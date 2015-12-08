allfiles = 1:50;
readdir = 'O:\Michael\2015\PSII supercomplex Pengqi\RUN III 20151001\experiment light treatment to rdouble reduce Qa\';
writedir = readdir;
APD = 1; % "1" is MPD, "0" is SPCM AQRH 16
filename2 = 'trace';
timeres = 0.1;
x=0.000:.004:15;
number='';

histdecayX=x(1:length(x)-1).';

parfor tr=allfiles
    
        complex=tr;
    
    [trace,delaytimes] = read_pt3_v4(timeres,readdir,[filename2 int2str(tr) number '.pt3'],APD);    
    
    histdecayY=hist(delaytimes(delaytimes>0),x); 
    histdecayY=histdecayY(1:length(histdecayY)-1).';
    histdecay=[histdecayX histdecayY];
        
    dlmwrite(fullfile(writedir,['trace' int2str(tr)]),trace,' ');
    dlmwrite(fullfile(writedir,['decaytrace' int2str(tr)]),histdecay,' ');
    tr
end


