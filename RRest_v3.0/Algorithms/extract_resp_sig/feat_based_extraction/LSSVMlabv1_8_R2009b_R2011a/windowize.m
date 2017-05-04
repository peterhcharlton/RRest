function w = windowize(A,window_array)
% Re-arrange the data points into a Hankel matrix for (N)AR time-series modeling
% 
% >> w = windowize(A, window)
%
% Use windowize function to make a nonlinear AR predictor with a
% nonlinear regressor. The last elements of the resulting matrix
% will contain the future values of the time-series, the others
% will contain the past inputs. window is the relative index of
% data points in matrix A, that are selected to make a window. Each
% window is put in a row of matrix W. The matrix W contains as many
% rows as there are different windows selected in A. 
% 
% Schematically, this becomes
% 
% >> A = [a1 a2 a3;
%          b1 b2 b3;
%          c1 c2 c3;
%          d1 d2 d3;  
%          e1 e2 e3;
%          f1 f2 f3;
%          g1 g2 g3];
% 
% >> W = windowize(A, [1 2 3])
% 
%    W = 
%      a1 a2 a3  b1 b2 b3  c1 c2 c3
%      b1 b2 b3  c1 c2 c3  d1 d2 d3  
%      c1 c2 c3  d1 d2 d3  e1 e2 e3  
%      d1 d2 d3  e1 e2 e3  f1 f2 f3 
%      e1 e2 e3  f1 f2 f3  g1 g2 g3
% 
% The function windowizeNARX converts the time-series and his
% exogeneous variables into a block hankel format useful for
% training a nonlinear function approximation as a nonlinear ARX
% model.  
% 
% Full syntax
%     (The length of window is denoted by w.)
% 
% >> Xw = windowize(X, window)
% 
%       Outputs    
%         Xw : (N-w+1) x w matrix of the sequences of windows over X
%       Inputs    
%         X  : N x 1 vector with data points
%         w  : w x 1 vector with the relative indices of one window
% 
%    
% see also:
%   windowizeNARX, predict, trainlssvm, simlssvm


% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab



l = max(window_array);
w = zeros(size(A,1)-l+1,length(window_array)*size(A,2));
for i=1:size(A,1)-l+1,
  for j = 1:length(window_array),
    w(i,(j-1)*size(A,2)+1:j*size(A,2)) = A(i-1+window_array(j),:);
  end
end
