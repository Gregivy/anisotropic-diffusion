function [ A,b,new_kappa ] = getAb( I,dt,kappa )

[n,m] = size(I);
A = zeros(n*m,5);
B = zeros(n+2,m+2);
b = zeros(1,n*m);

% tor like edges
B(1,2:m+1) = I(n,1:m);
B(1,1) = I(1,1);
B(1,m+2) = I(1,m);
B(n+2,2:m+1) = I(1,1:m);
B(n+2,1) = I(n,1);
B(n+2,m+2) = I(n,m);
B(2:n+1,1) = I(1:n,m);
B(2:n+1,m+2) = I(1:n,1);
B(2:n+1,2:m+1) = I;

g = @(x) 1/(1+(x/kappa)^2);

new_kappa = zeros(n*m,1);
k = 1;

for i=2:n+1
    for j=2:m+1
        fi = zeros(1,5);

        Qs = abs([B(i-1,j),B(i,j-1),B(i,j+1),B(i+1,j)] - B(i,j));

        minQ = min(Qs);
        maxQ = max(Qs);

        fi(1)=-dt*g(B(i-1,j)-B(i,j));
        fi(2)=-dt*g(B(i,j-1)-B(i,j));
        fi(3)=1+dt*(g(B(i-1,j)-B(i,j))+g(B(i,j-1)-B(i,j))+g(B(i+1,j)-B(i,j))+g(B(i,j+1)-B(i,j)));
        fi(4)=-dt*g(B(i,j+1)-B(i,j));
        fi(5)=-dt*g(B(i+1,j)-B(i,j));

        b((i-2)*m+j-1)=b((i-2)*m+j-1)+B(i,j);

        A(k,:) = fi;

        P = g(maxQ)/g(minQ);
        D = g(minQ)-g(maxQ);

        Npd = min(P,D);

        new_kappa(k) = maxQ*sqrt(Npd/(1-Npd));

        k=k+1;

    end
end

new_kappa = median(new_kappa);

end

