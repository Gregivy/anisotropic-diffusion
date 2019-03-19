function [ b ] = multAx( A,x,n,m )

k=1;
b = zeros(1,n*m);
for i=1:n
    for j=1:m
        s=0;
        %i-1,j
        if (i-1>=1 && i-1<=n && j>=1 && j<=m)
            s=s+x((i-2)*m+j)*A(k,1);
        else
            s=s+x((n-1)*m+j)*A(k,1);
        end
        %i,j-1
        if (i>=1 && i<=n && j-1>=1 && j-1<=m)
            s=s+x((i-1)*m+j-1)*A(k,2);
        else
            s=s+x((i-1)*m+m)*A(k,2);
        end
        %i,j
        s=s+x((i-1)*m+j)*A(k,3);
        %i,j+1
        if (i>=1 && i<=n && j+1>=1 && j+1<=m)
            s=s+x((i-1)*m+j+1)*A(k,4);
        else
            s=s+x((i-1)*m+1)*A(k,4);
        end
        %i+1,j
        if (i+1>=1 && i+1<=n && j>=1 && j<=m)
            s=s+x(i*m+j)*A(k,5);
        else
            s=s+x(j)*A(k,5);
        end
        b(k)=s;
        k=k+1;
    end
    
end

end

