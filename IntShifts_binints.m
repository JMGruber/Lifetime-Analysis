function newtrace = IntShifts_binints(trace,intbin)

binnedtrace=ones(1,floor(length(trace(:,1))/intbin));
traceint=trace(:,2);
for i=1:floor(length(trace(:,1))/intbin)
    binnedtrace(i)=sum(traceint((i-1)*intbin+1:i*intbin));
end
newtrace=[intbin*trace(1,1)*(1:length(binnedtrace))' binnedtrace'];