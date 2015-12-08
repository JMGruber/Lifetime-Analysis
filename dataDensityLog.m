function [ dmap ] = dataDensityLog( x, y, height, deltax, limits, fudge )
%DATADENSITY Get a data density image of data 
%   x, y - two vectors of equal length giving scatterplot x, y co-ords
%   height - y-dimension of the data density plot, in pixels 
%   deltax - x-scale of the data density plot, in pixels 
%           usually [0.01*(exp(0:log(1.5):9)) 120];
%   limits - [xmin xmax ymin ymax] - defaults to data max/min
%   fudge - the amount of smear, defaults to size of pixel diagonal
%
% By Malcolm McLean. Modified by TPJK.
%

    if(nargin < 4)
        deltax = [0.01*(exp(0:log(1.5):9)) 120];
    end
    
    limits(1) = 0.01;
    limits(2) = 100;
    limits(3) = min(y);
    limits(4) = max(y);
    
    deltay = (limits(4) - limits(3)) / height;
    
    if(nargin < 6)
        fudge = deltay/2;
    end

    width=length(deltax);
    dmap = zeros(height, width);
    for ii = 0: height - 1
        yi = limits(3) + ii * deltay + deltay/2;
        for jj = 1 : width
            xi = deltax(jj); %+ deltax/2;
            dd = 0;
            for kk = 1: length(x)
                dist2 = (log(x(kk)) - log(xi))^2 + (y(kk) - yi)^2;
                dd = dd + 1/(dist2 + fudge); 
            end
            dmap(ii+1,jj) = dd;
        end
    end
            

end