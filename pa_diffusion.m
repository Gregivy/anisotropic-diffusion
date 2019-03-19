function [ I ] = pa_diffusion( varargin )

% #1. pa_diffusion(I,dt,kappa) - uses article CG
% #2. pa_diffusion(I,dt,kappa,tol) - uses bult in CG with tolerance

if (nargin==3)
    I = varargin{1};
    dt = varargin{2};
    kappa = varargin{3};
elseif (nargin==4)
    I = varargin{1};
    dt = varargin{2};
    kappa = varargin{3};
    tol = varargin{4};
else
    error('pa_diffusion accepts up to 4 input arguments!')
end

[n,m] = size(I);

x = reshape(I',1,n*m);

if(nargin==4)
    x = x';
end

[A,b,kappa] = getAb(I,dt,kappa);
prevIm = x+2;

while(max(abs(x-prevIm))>1)

    prevIm = x;
    
    if (nargin==4)
        % #1. https://www.mathworks.com/help/matlab/ref/pcg.html
        matrixMult = @(v) multAx(A,v,n,m)';
        x = pcg(matrixMult,b',tol,[],[],[],x);
    else
        % #2. article proposal
        r = b - multAx(A,x,n,m);
        z = r;

        prevX = 10^10;

        while (max(abs(x-prevX))>1 && sum(r~=0))
            prevX = x;
            az = multAx(A,z,n,m);
            dotr = dot(r,r);
            alpha = dotr/dot(az,z);
            x = x + alpha*z;
            r = r -alpha*az;
            betta = dot(r,r)/dotr;
            z = r + betta*z;
        end
    end
    

    I = reshape(x,m,n)';
    [A,b,kappa] = getAb(I,dt,kappa);

end

end