function model = plotlssvm(model,ab,grain, princdim)
% Plot the LS-SVM results in the environment of the training data
% 
% >> plotlssvm({X,Y,type,gam, sig2, kernel})
% >> plotlssvm({X,Y,type,gam, sig2, kernel}, {alpha,b})
% >> model = plotlssvm(model)
% 
% The first argument specifies the LS-SVM. The latter specifies the
% results of the training if already known. Otherwise, the training
% algorithm is first called. One can specify the precision of the
% plot by specifying the grain of the grid. By default this value
% is 50. The dimensions (seldims) of the input data to display can
% be selected as an optional argument in case of higher dimensional
% inputs (> 2). A grid will be taken over this dimension, while the
% other inputs remain constant (0).
%  
%
% Full syntax
% 
%     1. Using the functional interface:
% 
% >> plotlssvm({X,Y,type,gam,sig2,kernel,preprocess}, {alpha,b})
% >> plotlssvm({X,Y,type,gam,sig2,kernel,preprocess}, {alpha,b}, grain)
% >> plotlssvm({X,Y,type,gam,sig2,kernel,preprocess}, {alpha,b}, grain, seldims)
% >> plotlssvm({X,Y,type,gam,sig2,kernel,preprocess})
% >> plotlssvm({X,Y,type,gam,sig2,kernel,preprocess}, [],      , grain)
% >> plotlssvm({X,Y,type,gam,sig2,kernel,preprocess}, [],      , grain, seldims)
% 
%       Inputs    
%         X             : N x d matrix with the inputs of the training data
%         Y             : N x 1 vector with the outputs of the training data
%         type          : 'function estimation' ('f') or 'classifier' ('c')
%         gam           : Regularization parameter
%         sig2          : Kernel parameter (bandwidth in the case of the 'RBF_kernel')
%         kernel(*)     : Kernel type (by default 'RBF_kernel')
%         preprocess(*) : 'preprocess'(*) or 'original'
%         alpha(*)      : support values obtained from training
%         b(*)          : Bias term obtained from training
%         grain(*)      : The grain of the grid evaluated to compose the surface (by default 50)
%         seldims(*)    : The principal inputs one wants to span a grid (by default [1 2])
% 
%
%     2. Using the object oriented interface:
% 
% >> model = plotlssvm(model)
% >> model = plotlssvm(model, [], grain)
% >> model = plotlssvm(model, [], grain, seldims)
% 
%       Outputs    
%         model(*)   : Trained object oriented representation of the LS-SVM model
%       Inputs    
%         model      : Object oriented representation of the LS-SVM model
%         grain(*)   : The grain of the grid evaluated to compose the surface (by default 50)
%         seldims(*) : The principal inputs one wants to span a grid (by default [1 2])
% 
% See also:
%   trainlssvm, simlssvm.

% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab

fprintf('Start Plotting...')

%
% initiating the model...
%
if iscell(model),   
  model = initlssvm(model{:});  
  eval('model.alpha = ab{1}; model.b = ab{2};model.status = ''trained'';','model=trainlssvm(model);');
end



%figure;
clf

model = trainlssvm(model);
% reconstruct the original support vectors ...
[osvX,osvY] = postlssvm(model,model.xtrain(:,1:model.x_dim),model.ytrain(:,1:model.y_dim));


%
% define the principal dimensions one plots
%
if (model.x_dim>2) 
  % plotted principal dimensions
  eval('princdim; restdim = setdiff(1:model.x_dim,princdim);','princdim=[1 2 3];');
elseif (model.x_dim==2),
  princdim = [1 2]; restdim = []; 
else
  princdim = [1]; restdim = []; 
end

if max(princdim)>model.x_dim, 
  error('Given dimensions exceed input dimensions...');
end


% classification (x_dim=2, y_dim=1:...) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if model.type(1)=='c', % 'classification' 

  %
  % precision of plot
  %
  eval('grain;','grain = 150;');
  
  if model.x_dim>=2, 
  %%%%%%%%%%%%%%%%%%
   
    % Determine plot limits 
    xmin1=min(osvX(:,princdim(1))); if xmin1<0, xmin1=1.05*xmin1; else xmin1 = 0.98*xmin1; end
    xmax1=max(osvX(:,princdim(1))); if xmax1>0, xmax1=1.05*xmax1; else xmax1 = 0.98*xmax1; end
    xmin2=min(osvX(:,princdim(2))); if xmin2<0, xmin2=1.05*xmin2; else xmin2 = 0.98*xmin2; end
    xmax2=max(osvX(:,princdim(2))); if xmax2>0, xmax2=1.05*xmax2; else xmax2 = 0.98*xmax2; end
    xrange1 = xmin1:(xmax1-xmin1)/grain:xmax1;
    xrange2 = xmin2:(xmax2-xmin2)/grain:xmax2;
    [XX,YY] = meshgrid(xrange1,xrange2);
    Xt = [reshape(XX,numel(XX),1) reshape(YY,numel(YY),1)];
    xsteps = length(xrange1);
    ysteps = length(xrange2);
    
    
    
    %
    % simulate the points
    %
    restdim = setdiff(1:model.x_dim, princdim);
    rest = zeros(size(Xt,1),model.x_dim-2);
    Xt = [Xt rest];
    [ZZ,~,model] = simlssvm(model,Xt(:,[princdim restdim]));
    if min(ZZ)==max(ZZ), warning('Simulation over the input space results in only one class...'); end
    
    %
    % for plotting, the categorical format is required
    %
    if ~strcmpi(model.codetype,'none'),
      if size(model.codebook1,1)~=1,
	eval('[ZZ,codebook_cat] = code(ZZ,''code_cat'',model.codebook2,model.code_distfct);',...
	     '[ZZ,codebook_cat] = code(ZZ,''code_cat'',model.codebook2);');
      else
	codebook_cat = model.codebook1;
      end
      eval('osvY = code(osvY, codebook_cat,{}, model.codebook2, model.codedist_fct, model.codedist_args);',...
	   'osvY = code(osvY, codebook_cat,{}, model.codebook2);');    
      if max(max(ZZ))==-inf, 
	error('bad coding scheme, no classes found after training');
      end

    else
      
      if model.y_dim>1,
	warning(['only first dimension is plotted, for multiclass' ...
		 ' classification use categorical representation, ev.'...
		 ' combined with a coding technique.']);
      end
      osvY = osvY(:,1);
      ZZ = ZZ(:,1);
      sosvY = sort(osvY);
      codebook_cat = sosvY([1;find(sosvY(2:end)~=sosvY(1:end-1))+1])';
    end

    
    % contour plot
    colormap cool;
    map = colormap;
    %cindex = [min(codebook1)+.1 codebook1 max(codebook1)-.1];
    ZZd = reshape(ZZ(:,1),size(XX,1),size(XX,2));
    eval('[C,h]=contourf(XX,YY,ZZd);','warning(''no surface plot feasable'');'); 
    hold on;
    eval('clabel(C,h,codebook_cat);',' ');
    
    
    %
    % plotting the datapoints
    %
    legstr{1} = 'Classifier';
    markers = {'*','s','+','o','x','d','v','p','h'};
    for c=1:length(codebook_cat),
      s = find(osvY(:,1)==codebook_cat(c));

      plot(osvX(s,princdim(1)),osvX(s,princdim(2)) ,[markers{1+mod(c-1,9)} 'k']);
      legstr{c+1} = ['class ' num2str(c)];
    end
    eval('legend(legstr{1:end});');
    
    
    % arrange axis
    xlabel(['X_{' num2str(princdim(1)) '}']);
    ylabel(['X_{' num2str(princdim(2)) '}']);
    title(['LS-SVM_{\gamma=' num2str(model.gam(1)) ',\sigma^2=' num2str(model.kernel_pars(1:min(1,length(model.kernel_pars)))) '}^{' model.kernel_type(1:3) '}, with ' num2str(length(codebook_cat)) ' different classes']);
    axis([xmin1 xmax1 xmin2 xmax2]);  
    hold off;    
  else    
    error('cannot display this dimension..');
  end
  
  
  
% function estimation (x_dim=1,2; y_dim=1)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif model.type(1)=='f',
  eval('grain;','grain = 200;');
  
  
  % Determine plot limits   
  xmin1=min(osvX(:,princdim(1))); 
  xmax1=max(osvX(:,princdim(1))); 
    
  if model.x_dim>=2 & length(princdim)==2,
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Determine plot limits 
    xmin2=min(osvX(:,princdim(2))); if xmin2<0, xmin2=1.05*xmin2; else xmin2 = 0.975*xmin2; end
    xmax2=max(osvX(:,princdim(2))); if xmax2>0, xmax2=1.05*xmax2; else xmax2 = 0.975*xmax2; end
    range1 = (xmin1:(xmax1-xmin1)/grain:xmax1)';
    range2 = (xmin2:(xmax2-xmin2)/grain:xmax2)';
    
    rest = zeros(size(range1,1),model.x_dim-2);
    for i=1:length(range2),
      Xt = [range1 ones(size(range1,1)).*range2(i) rest];
      [r,~,model]  = simlssvm(model, Xt(:,[princdim,restdim]));
      z(i,:)=r';
    end    
    surf(range1, range2,z);
    hold on;
    plot3(osvX(model.selector,princdim(1)),osvX(model.selector,princdim(2)), osvY(model.selector,1),'k*');
    shading interp;

    xlabel(['X_' num2str(princdim(1))]);
    ylabel(['X_' num2str(princdim(2))]);
    zlabel('Y');
    title([' function estimation using LS-SVM_{\gamma=' num2str(model.gam(1)) ',\sigma^2=' num2str(model.kernel_pars(1:min(1,length(model.kernel_pars)))) '}^{' model.kernel_type(1:3) '} ']);
    view(-30,50);
    hold off;
  elseif and(model.x_dim==1,model.y_dim==1) || length(princdim)==1,

    range1 = (xmin1:(xmax1-xmin1)/grain:xmax1)'; 
    rest = zeros(size(range1,1),model.x_dim-1);
    grid = [range1 rest];
    [z,~,model]  = simlssvm(model,grid(:,[princdim(1) restdim]) ); 
    plot(range1,z,'r-');
    hold on;
    plot(osvX(model.selector,princdim(1)),osvY(model.selector,1),'b.');

    xlabel('X');
    ylabel('Y');
    if model.kernel_type(1)~='p'
        title([' function estimation using  LS-SVM_{\gamma=' num2str(model.gam(1)) ',\sigma^2=' num2str(model.kernel_pars(1:min(1,length(model.kernel_pars)))) '}^{' model.kernel_type(1:3) '}']);
    else
        title([' function estimation using  LS-SVM_{\gamma=' num2str(model.gam(1)) ',t=' num2str(model.kernel_pars(1)) ',degree=' num2str(model.kernel_pars(2)) '}^{' model.kernel_type(1:3) '}']);
    end
    %eval('title(['' function estimation using  LS-SVM_{\gamma='' num2str(model.gam(1)) '',\sigma^2='' num2str(model.kernel_pars) ''}^{'' kerneltype ''} datapoints (black *), and estimation  (blue line)'']);',' title(''function approximation using LS-SVM'')');

    hold off;

  else
    Yh = simlssvm(model,osvX);
    plot(Yh,'r-');
    hold all; 
    plot(osvY,'.');
    xlabel('time');
    ylabel('Y');
    title([' function estimation using '...
	   ' LS-SVM_{\gamma=' num2str(model.gam(1)) ...
	   ',\sigma^2=' num2str(model.kernel_pars(1:min(1,length(model.kernel_pars)))) '}^{' model.kernel_type(1:3) '}'...
	   ' datapoints (blue .), and estimation  (red line)']);
    hold off

  end
  
else
  
  
end


fprintf('finished\n');
