function [wX, wY, xdim, ydim, n] = windowizeNARX(X,Y,x_delays,y_delays, steps)
% Re-arrange the data points into a block Hankel matrix for (N)ARX time-series modeling
% 
% >> [Xw,Yw] = windowizeNARX(X,Y,xdelays, ydelays, steps)
% 
%  Rearrange the points of X and Y in a regressor matrix of the
%  past inputs and outputs (Xw) and the future outputs (Yw). 
% 
% Full syntax
%
% >> [Xw, Yw, xdim, ydim, n] = windowizeNARX(X, Y, xdelays, ydelays, steps)
% 
%       Outputs    
%         Xw Matrix of the data used for input including the delays
%         Yw Matrix of the data used for output including the next steps
%         xdim(*) Number of dimensions in new input
%         ydim(*) Number of dimensions in new output
%         n(*) Number of new data points
%       Inputs    
%         X       : N x m vector with input data points
%         Y       : N x d vector with output data points
%         xdelays : Number of lags of X in new input
%         ydelays : Number of lags of Y in new input
%         steps(*): Number of future steps of Y in new output (by default 1)
%
% See also:
%   windowize, predict, trainlssvm, simlssvm

% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab


m=max(x_delays,y_delays);

eval('steps;','steps = 1;');
if steps == 0, 
  n = size(X,1)-m;
else
  n = size(X,1)-m -steps+1;
end


wX = zeros(n,size(X,2)*(x_delays+1)+size(Y,2)*(y_delays));
wY = zeros(n,size(Y,2)*steps);
xdim = size(wX,2);
ydim = size(wY,2);

hdx  = (x_delays+1)*size(X,2);


for t=1:n,
  for i=1:x_delays+1,
    wX(t,1+((i-1)*size(X,2):i*size(X,2)-1)) = X(t+m-x_delays+i-1,:);
  end
  
  for i=1:y_delays,
    wX(t,hdx + (i-1)*size(Y,2) + (1:size(Y,2))) = Y(t+m-y_delays+i-1,:);   
  end

  for i=1:steps,
    wY(t,i:i+size(Y,2)-1) = Y(t+m+i-1,:);   
  end
  
end  
  
  