function [model,b] = robustlssvm(model,ab,X,Y)
% Robust training in the case of non-Gaussian noise or outliers
%(only possible with the object oriented interface)
%
% >> model     = robustlssvm(model)
%
% Robustness towards outliers can be achieved by reducing the
% influence of support values corresponding to large errors.
%
%
% Full syntax
%
%     1. Using the object oriented interface:
%
% >> model = robustlssvm(model)
%
%       Outputs
%         model : Robustly trained object oriented representation of the LS-SVM model
%       Inputs
%         model : Object oriented representation of the LS-SVM model
%
% See also:
%   trainlssvm, tunelssvm, crossvalidate


% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab



if iscell(model),
    func = 1;
    model = initlssvm(model{:});
else
    func = 0;
end


if model.type(1)~='f',
    error('Robustly weighted least squares only implemented for regression case...');
end



if nargin>1,
    if iscell(ab) && ~isempty(ab),
        model.alpha = ab{1};
        model.b = ab{2};
        model.status = 'trained';
        if nargin>=4,
            model = trainlssvm(model,X,Y);
        end
    else
        model = trainlssvm(model,ab,X);
    end
else
    model = trainlssvm(model);
end


% model errors
ek = model.alpha./model.gam';
g = model.gam;
%
% robust estimation of the variance
%
eval('delta=model.delta;','delta=[];')
for j=1:500
    vare = 1.483*median(abs((ek)-median(ek)));
    alphaold = model.alpha;
    %
    % robust re-estimation of the alpha's and the b
    %
    cases = reshape((ek./vare),1,model.nb_data);
    W = weightingscheme(cases,model.weights,delta);
    W = g*W;
    
    model = changelssvm(model,'gam',W);
    %     model = changelssvm(model,'implementation','MATLAB');
    model = trainlssvm(model);
    ek = model.alpha./model.gam';
    
    if norm(abs(alphaold-model.alpha),'fro')<=1e-4,
        fprintf('Converged after %.0f iteration(s)', j);
        if func && nargout~=1,
            b = model.b;
            model = model.alpha;
        end
        return
    end
    model.status = 'changed';
end

