allfiles = 20:30;
readdir = 'C:\CD\Data WorkingOn\pH 4\';
writedir = 'C:\CD\Data WorkingOn\pH 4\';
filename2 = 'trace';
channel = 2; % 1 for new APD, 2 for old APD
timeres = 0.01;
x=.000:.004:13.15;

histdecayX=x(1:length(x)-1).';

for tr=allfiles
    
        complex=tr;
    
    [trace,delaytimes] = read_pt3_v4_channel(timeres,readdir,[filename2 int2str(tr) '.pt3'],channel);    
    
    histdecayY=hist(delaytimes(delaytimes>0),x); 
    histdecayY=histdecayY(1:length(histdecayY)-1).';
    histdecay=[histdecayX histdecayY];
        
    dlmwrite(fullfile(writedir,['trace' int2str(tr)]),trace,' ');
    dlmwrite(fullfile(writedir,['Decaytrace' int2str(tr)]),histdecay,' ');
	 fprintf('%s%6.2f\n','trace',tr);
end


