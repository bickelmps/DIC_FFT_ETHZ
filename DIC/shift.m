function [B]=shift(A,nr,nc)
% function [B]=shift[A,nx,ny]
%
% Shifts the matrix A of a certain number of rows (nr) and colums (nc)
% positive and negative values

[M,N]=size(A);

B=zeros(M,N);

if nr>=0 & nc>=0
    
    B(nr+1:M,nc+1:N)=A(1:M-nr,1:N-nc);

elseif nr<0 & nc<0
    
    B(1:M+nr,1:N+nc)=A(1-nr:M,1-nc:N);
    
elseif nr>=0 & nc<0
    
    
    B(nr+1:M,1:N+nc)=A(1:M-nr,1-nc:N);
    
    
elseif nr<0 & nc>=0
   
    B(1:M+nr,nc+1:N)=A(1-nr:M,1:N-nc);
        
    
end
    


return