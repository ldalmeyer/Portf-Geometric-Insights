function [b] =compound_returns(x)
a=1;
    for i=1:length(x);
        x_i=x(i);
       
        a=(a*(1+x_i));
       
    end
    
    b=a-1;
end