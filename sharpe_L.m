function [sr] = sharpe_L(x) 
   ncols = size(x);
   if istable(x)
    x = table2array(x);
   end
   er=[];
   sr=[];
   for i=1:ncols(2)
    er = x(:,i);
    sr(i) = (mean(er)/std(er));
   end
end

