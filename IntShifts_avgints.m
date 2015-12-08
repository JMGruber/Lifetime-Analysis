function newtrace = IntShifts_avgints(trace,intbin)

% This function averages the number "intbin" of consecutive intensities

binnedtrace=ones(1,floor(length(trace(:,1))/intbin));
traceint=trace(:,2);
for i=1:floor(length(trace(:,1))/intbin)
    b=traceint((i-1)*intbin+1:i*intbin);
    binnedtrace(i)=sum(b)./length(b);
end
newtrace=[intbin*trace(1,1)*(1:length(binnedtrace))' binnedtrace'];