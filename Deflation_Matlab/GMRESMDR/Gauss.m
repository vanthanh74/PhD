function [X] = Gauss(A)
m=size(A,1);
n=size(A,2);
i=1;
j=1;
while( i<= m && j <= n )    
        if ( A(i,j) == 0)                  
            dem=0;
            for k=i+1:m             
                  if(A(k,j)~=0)                    
                      dem = dem + 1;                                           
                  end
            end
             if(dem == 0)
                j=j+1;
             else 
                 for k=i+1:m             
                      if(A(k,j)~=0)                     
                             Temp = A(i,:);
                             A(i,:)=A(k,:);
                             A(k,:)=Temp; 
                             break;
                      end
                 end  
             end
        else
            for k=i+1:m
                A(k,:)=A(k,:)-A(k,j)/A(i,j)*A(i,:);
            end
            i=i+1;
            j=j+1;
        end  
end
X=A;