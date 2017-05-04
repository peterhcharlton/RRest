function [Xm,Xval,itr,fig] = linesearch(fun,startvalues,funargs,varargin)
% Global optimization by exhaustive search over the parameter space   
%
% >> Xopt = gridsearch(fun, startvalues)
%
% The most simple algorithm to determine the minimum of a cost
% function with possibly multiple optima is to evaluate a grid over
% the parameter space and to pick the minimum. This procedure
% iteratively zooms to the candidate optimum.
% The startvalues determine the limits of the grid over parameter
% space. 
%
%
% Full syntax
%
% >> [Xopt, Yopt, Evaluations, fig] = linesearch(fun, startvalues, funargs, option1,value1,...)
% 
%       Outputs    
%         Xopt           : Optimal parameter set
%         Yopt           : Criterion evaluated at Xopt
%         Evaluations    : Used number of function evaluations
%         fig            : handle to used figure
%       Inputs    
%         fun Function   : implementing the cost criterion
%         startvalues    : 2*d matrix with starting values of the optimization routine
%         funargs(*)     : Cell with optional extra function arguments of fun
%         option (*)     : The name of the option one wants to change
%         value  (*)     : The new value of the option one wants to change
%
% The different options and their meanings are:
%
%         Nofigure      : 'figure'(*) or 'nofigure'
%         MaxFunEvals   : Maximum number of function evaluations (default: 20)
%         GridReduction : grid reduction parameter (e.g. '1.5':
%                           small reduction; `10': heavy reduction; default '2')
%         TolFun        : Minimal toleration of improvement on function value (default: 0.01) 
%         TolX          : Minimal toleration of improvement on X value (default: 0.01)
%         Grain         : number of evaluations per iteration (default: 10)
% 
% see also:
%   gridsearch, tunelssvm, crossvalidate, fminunc
  
  
% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab


dim = prod(size(startvalues));
if dim~=2,
  error('optimization only possible for 1-dimensional problems...');
end


%
% defaults
%
nofigure='figure';
eval('funargs;','funargs={};');
maxFunEvals = 20;
TolFun = .01;
TolX = .01;
grain = 10;
zoomfactor=2;


%
% extra input arguments
%
for t=1:2:length(varargin),
  if t+1>length(varargin),
    warning('extra arguments should occur in pairs (''name'',value)');
  else  
    if strcmpi(varargin{t},'nofigure'),        nofigure=varargin{t+1};
    elseif strcmpi(varargin{t},'maxFunEvals'), maxFunEvals=varargin{t+1}
    elseif strcmpi(varargin{t},'zoomfactor'),  zoomfactor=varargin{t+1}
    elseif strcmpi(varargin{t},'TolFun'),      TolFun=varargin{t+1}
    elseif strcmpi(varargin{t},'TolX'),        TolX=varargin{t+1}
    elseif strcmpi(varargin{t},'grain'),       grain=varargin{t+1}; 
    else  warning(['option ' varargin{t} ' unknown']);
    end
  end
end

itr =maxFunEvals;


%
% initiate grid
%
graina = -.5:1/(grain-1):.5;
grid   = [min(startvalues) max(startvalues)];
center = mean(grid);

if nofigure(1) =='f',fig = figure; hold on; end


itr = 0;
Xm_old = inf;
Xval_old = inf;
Xm = -inf;
Xval = -inf;
while itr<maxFunEvals & norm(Xm-Xm_old)>TolX & norm(Xval-Xval_old)>TolFun,
  Xm_old   = Xm;
  Xval_old = Xval;
  
  xtrma = [min(startvalues) max(startvalues)];
  xline = xtrma(1):(xtrma(2)-xtrma(1))/(grain-1):xtrma(2);
  for i = 1:length(xline),
    cost(i) = feval(fun, xline(i), funargs{:});
    itr = itr+1;
    if nofigure(1) =='f',
      plot(xline(i),cost(i),'dk');drawnow
    end
  end

  [sc, si] = sort(cost);
  Xm   = xline(si(1)); 
  Xval = sc(1); 
  selected = si(1:ceil(length(si)/zoomfactor));
  startvalues = [min(xline(selected)) max(xline(selected))];
  
end
  
  


