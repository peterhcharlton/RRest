function [Xm,Xval,evals,fig] = gridsearch(fun,startvalues,funargs,varargin)
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
% This optimisation function can be customized by passing extra
% options and the corresponding value
%
% >> [Xopt, Yopt, Evaluations, fig] = gridsearch(fun, startvalues, funargs, option1,value1,...)
%
% the possible options and their default values are: 
%  'nofigure'   ='figure';
%  'maxFunEvals'= 190;
%  'TolFun'     = .0001;
%  'TolX'       = .0001;
%  'grain'      = 10;
%  'zoomfactor' =5;
%
%
% An example is given:
%
% >> fun = inline('1-exp(-norm([X(1) X(2)]))','X');
% >> gridsearch(fun,[-3 3; 3 -3])
%
% or 
% >> gridsearch(fun,[-3 3; 3 -3],{},'nofigure','nofigure','MaxFunEvals',1000)
%
%
% Full syntax
%
% >> [Xopt, Yopt, Evaluations, fig] = gridsearch(fun, startvalues, funargs, option1,value1,...)
% 
%       Outputs    
%         Xopt           : Optimal parameter set
%         Yopt           : Criterion evaluated at Xopt
%         Evaluations    : Used number of cost function evaluations
%         fig            : handle to used figure
%       Inputs    
%         CostFunction   : implementing the cost criterion
%         Startvalues    : 2*d matrix with starting values of the optimization routine
%         Funargs(*)     : Cell with optional extra function arguments of fun
%         option (*)     : The name of the option one wants to change
%         value  (*)     : The new value of the option one wants to change
%
% The different options and their meanings are:
%
%         Nofigure      : 'figure'(*) or 'nofigure'
%         MaxFunEvals   : Maximum number of function evaluations (default: 100)
%         GridReduction : grid reduction parameter (e.g. '2':
%                           small reduction; `10': heavy reduction; default '5')
%         TolFun        : Minimal toleration of improvement on function value (default: 0.0001) 
%         TolX          : Minimal toleration of improvement on X value (default: 0.0001)
%         Grain         : Square root number of function evaluations in one grid (default: 8)
% 
% 
% see also:
%   tunelssvm, crossvalidate, fminunc
  
  
% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab


dim = size(startvalues,2);
if dim~=2,
  error('optimization only possible for 2-dimensional problems...');
end


%
% defaults
%
nofigure='nofigure';
eval('funargs;','funargs={};');
maxFunEvals = 70;
TolFun = .0001;
TolX = .0001;
grain = 7;
zoomfactor=5;


%
% extra input arguments
%
for t=1:2:length(varargin),
  if t+1>length(varargin),
    warning('extra arguments should occur in pairs (''name'',value)');
  else  
    if strcmpi(varargin{t},'nofigure'),        nofigure=varargin{t+1};
    elseif strcmpi(varargin{t},'maxFunEvals'), maxFunEvals=varargin{t+1};
    elseif strcmpi(varargin{t},'zoomfactor'),  zoomfactor=varargin{t+1};
    elseif strcmpi(varargin{t},'TolFun'),      TolFun=varargin{t+1};
    elseif strcmpi(varargin{t},'TolX'),        TolX=varargin{t+1};
    elseif strcmpi(varargin{t},'grain'),       grain=varargin{t+1}; 
    else  warning(['option ' varargin{t} ' unknown']);
    end
  end
end

itr =maxFunEvals;
grain

%
% initiate figure
%
if nofigure(1) =='f', fig = figure;end
marx = {'o','x','+','*','s','d','v','^','<','>','p','h'};

%
% initiate output
%
disp('FIRST ITERATION:');


%
% initiate grid
%
dirs = [1 0; 0 1];
graina = -.5:1/(grain-1):.5;
for d=1:dim, 
  grid(d,:) = [min(startvalues(:,d)) max(startvalues(:,d)) 0];
end
center = mean(grid(:,1:2),2)';


for i=1:grain,
  for j = 1:grain,
    gridF((i-1)*grain+j,1:2) = center+dirs(1,:).*(grid(1,2)-grid(1,1)).*graina(i)+...
	dirs(2,:).*(grid(2,2)-grid(2,1)).*graina(j);
    gridF((i-1)*grain+j,3) = feval(fun, gridF((i-1)*grain+j,1:2),funargs{:});
    itr = itr-1;

    % figure...
    if nofigure(1) =='f',
      plot3(gridF(:,1),gridF(:,2),gridF(:,3),'.k');
      hold on;
      plot3(gridF((i-1)*grain+j,1),gridF((i-1)*grain+j,2),gridF((i-1)*grain+j,3),'bo','linewidth',5);
      title(['computing element [' num2str([i j]) '] of [' num2str([grain grain]) ']...']);hold off;
      drawnow; 
    end
  end
end


if nofigure(1) =='f',
  plot3(gridF(:,1),gridF(:,2),gridF(:,3),'.k');
  hold on; 
  contour(gridF(1:grain:end,1),gridF(1:grain,2),reshape(gridF(:,3),grain,grain),grain);  
  title('Optimization grid after first iteration');
  %view(-22,50); 
  view(0,90); 
  drawnow; 
end



zoom = range(gridF(:,1))*range(gridF(:,2));

% init optims
[optF(1),opti] = min(gridF(:,3),[],1);
optX(1,:) = gridF(opti,1:2);
dF = inf;
dX = inf;

disp(['   X=' num2str(optX(1,:)) ' ,F(X)=' num2str(optF(1)) ';']);


% indexes of X1 and X2
x1g = 1:grain:grain*grain;
x2g = 1:grain;

gridFold = [];
t=2;  
while and(itr>0, and(dF>TolFun, dX>TolX)), 
  disp(['ITERATION: ' num2str(t)]);

  %
  % re-new startvalues
  %
  center = optX(t-1,:);
  oldF = optF(t-1);
  
  
    
  % optimal directions: trimed mean variant
  % 1. select the alpha~1/zoomfactor least samples.
  [~,si] = sort(gridF(:,3));
  sis = si(1:ceil(grain*grain/zoomfactor));
  
  % 2. build a new grid based on this selected
  gridR = [gridF(sis,1)-center(1) gridF(sis,2)-center(2)];
  % optimal directions V: norm V==1
  [~,S,dirs] = svd(gridR);
  S = diag(S).^.5;
  %line([center(1) center(1)+.5*S(1)*dirs(1,1)],[center(2) center(2)+.5*S(1)*dirs(2,1)]);
  %line([center(1) center(1)+.5*S(2)*dirs(1,2)],[center(2) center(2)+.5*S(2)*dirs(2,2)]);

  
  %
  % re-initiate grid and evaluate
  %
  gridFold = [gridFold;gridF];
  for i=1:grain,
    for j = 1:grain,
      gridF((i-1)*grain+j,1:2) = center+dirs(1,:)*S(1)*graina(i)+dirs(2,:)*S(2)*graina(j); 
      gridF((i-1)*grain+j,3) = feval(fun, gridF((i-1)*grain+j,1:2),funargs{:});
      itr = itr-1;

      % figure...
      if nofigure(1) =='f',
	hold on
	plot3(gridF((i-1)*grain+j,1),gridF((i-1)*grain+j,2), gridF((i-1)*grain+j,3),['k' marx{t}]);
	title(['computing element [' num2str([i j]) '] of [' num2str([grain grain]) ']...']);hold off;
	drawnow; 
      end

    end
  end
  
  % init optims
  [optF(t),opti] = min(gridF(:,3),[],1);
  optX(t,:) = gridF(opti,1:2);
  if oldF<optF(t),
    % if worse: contract 
    optF(t) = mean([oldF;optF(t)]);
    optX(t,:) = mean([center;optX(t,:)],1);
  end
  dF = abs(oldF-optF(t));
  dX = norm(center-optX(t,:));  
  
  
  
  disp(['  dF=' num2str(dF) ', dX=' num2str(dX) ', X=' num2str(optX(t,:)) ...
        ' ,F(X)=' num2str(optF(t)) ';']);
  t=t+1;
  
        
        end

%
% output
%
Xm=optX(t-1,:);
Xval = optF(t-1);
evals = maxFunEvals-itr;