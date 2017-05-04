function W = weightingscheme(cases,wfct,varargin)

%
% INTERNAL FUNCTION
%

% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab

switch wfct
    
    case {'WHuber','whuber'}
        b = varargin{1};
        cases = abs(cases);
        temp = 1;
        W = temp.*(cases<b) + temp.*(cases>=b).*(b./cases);
        
    case {'WLogistic','wlogistic'}
        W = tanh(cases)./cases;
        
    case {'WHampel','whampel'}
       cases = abs(cases);
       % defaults for c1 and c2
       c1=2.5; c2=3; 
       %[c1,c2] = adaptweight(cases);
       dc = c2-c1;
       temp = 1; 
       W = temp.*(cases<=c1) + ...
           temp.*(cases<=c2 & cases>c1).*((c2-cases)./dc) + ...
           temp.*(cases>c2).*10e-8;
       
    case {'WMyriad','wmyriad'}
        %K = 0.5*iqr(cases);
        K = varargin{1};
        W = K^2./(K^2+cases.^2);  
end