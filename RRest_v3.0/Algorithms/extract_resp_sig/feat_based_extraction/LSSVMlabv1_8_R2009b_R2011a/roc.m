function [AREA,SE,RESULT_S,FPR_ROC,TPR_ROC,TNa,TPa,FNa,FPa]=roc(RESULT,CLASS,fig)
% Receiver Operating Characteristic (ROC) curve of a binary classifier
% 
% >> [area, se, deltab, oneMinusSpec, sens, TN, TP, FN, FP] = roc(Zt, Y)
% 
% The ROC curve shows the separation abilities of a binary
% classifier: by iteratively setting the possible classifier
% thresholds, the dataset is tested on misclassifications. As a
% result, a plot is shown where the various outcomes are
% described. If the plot has a surface of 1 on test data, a
% perfectly separating classifier is found (on that particular
% dataset), if the area equals 0.5, the classifier has no
% discriminative power at all. In general, this function can be
% called with the latent variables Zt and the corresponding class labels Yclass
% 
% >> Zt       = [-.7             Yclass = [-1
%                 .3                       -1
%                1.5                        1 
%                ...                       ..   
%                -.2]                       1]
% >> roc(Zt, Yclass)
% 
% For use in LS-SVMlab, a shorthand notation allows making the ROC
% curve on the training data. Implicit training and simulation of
% the latent values simplifies the call.
% 
% >> roc({X,Y,'classifier',gam,sig2,kernel})
% >> roc(model)
% 
%
% Full syntax
% 
%     1. Standard call (LS-SVMlab independent):
% 
% >> [area, se, deltab, oneMinusSpec, sens, TN, TP, FN, FP] = roc(Zt, Y)
% >> [area, se, deltab, oneMinusSpec, sens, TN, TP, FN, FP] = roc(Zt, Y, figure)
% 
%       Outputs    
%         area(*)   : Area under the ROC curve
%         se(*)     : Standard deviation of the residuals
%         deltab(*) : N x 1 different thresholds value
%         oneMinusSpec(*) : 1-Specificity of each threshold value
%         sens(*)   : Sensitivity for each threshold value
%         TN(*)     : Number of true negative predictions
%         TP(*)     : Number of true positive predictions
%         FN(*)     : Number of false negative predictions
%         FP(*)     : Number of false positive predictions
%       Inputs    
%         Zt        : N x 1 latent values of the predicted outputs
%         Y         : N x 1 of true class labels
%         figure(*) : 'figure'(*) or 'nofigure'
% 
%
%     2. Using the functional interface for the LS-SVMs:
% 
% >> [area, se, deltab, oneMinusSpec, sens, TN, TP, FN, FP] = roc({X,Y,'classifier',gam,sig2,kernel})
% >> [area, se, deltab, oneMinusSpec, sens, TN, TP, FN, FP] = roc({X,Y,'classifier',gam,sig2,kernel}, figure)
% 
%       Outputs    
%         area(*)   : Area under the ROC curve
%         se(*)     : Standard deviation of the residuals
%         deltab(*) : Different thresholds
%         oneMinusSpec(*) : 1-Specificity of each threshold value
%         sens(*)   : Sensibility for each threshold value
%         TN(*)     : Number of true negative predictions
%         TP(*)     : Number of true positive predictions
%         FN(*)     : Number of false negative predictions
%         FP(*)     : Number of false positive predictions
%       Inputs    
%         X             : N x d matrix with the inputs of the training data
%         Y             : N x 1 vector with the outputs of the training data
%         type          : 'function estimation' ('f') or 'classifier' ('c')
%         gam           : Regularization parameter
%         sig2          : Kernel parameter (bandwidth in the case of the 'RBF_kernel')
%         kernel(*)     : Kernel type (by default 'RBF_kernel')
%         preprocess(*) : 'preprocess'(*) or 'original'
%         figure(*)     : 'figure'(*) or 'nofigure'
% 
%
%     3. Using the object oriented interface for the LS-SVMs:
% 
% >> [area, se, deltab, oneMinusSpec, sens, TN, TP, FN, FP] = roc(model)
% >> [area, se, deltab, oneMinusSpec, sens, TN, TP, FN, FP] = roc(model, figure)
% 
%       Outputs    
%         area(*)   : Area under the ROC curve
%         se(*)     : Standard deviation of the residuals
%         deltab(*) : N x 1 vector with different thresholds
%         oneMinusSpec(*) 1-Specificity of each threshold value
%         sens(*)   : Sensibility for each threshold value
%         TN(*)     : Number of true negative predictions
%         TP(*)     : Number of true positive predictions
%         FN(*)     : Number of false negative predictions
%         FP(*)     : Number of false positive predictions
%       Inputs    
%         model     : Object oriented representation of the LS-SVM model
%         figure(*) : 'figure'(*) or 'nofigure'
% 
% See also:
%   deltablssvm, trainlssvm

% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @ http://www.esat.kuleuven.be/sista/lssvmlab

eval('fig;','fig=''figure'';');


%
% roc(model)
%
if iscell(RESULT), 
  RESULT = initlssvm(RESULT{:});
end
if isstruct(RESULT),
  model = RESULT;
  
  if model.type~='c',
    error(' ROC only possible for classification...');
  end

  RESULT = latentlssvm(model,postlssvm(model,model.xtrain));
  
  if size(RESULT,2)>1,
    warning(' ROC only possible for binary classification...');
  end  
  
  CLASS = codelssvm(model,model.ytrain);
  
end

if min(CLASS)~=-1 | max(CLASS)~=1,
    warning('Class labels need to be -1 or 1. Automatic conversion to -1 and 1 ...');
    s = find(CLASS==0); CLASS(s)=-1;  
end

if numel(CLASS)~=numel(RESULT),
  warning('Number of elements in Zt and Y must be equal;');
end



%
% roc(RESPONS, CLASS)
%

FI=find(isfinite(RESULT));
RESULT=(RESULT(FI));
CLASS=CLASS(FI);

NRSAM=size(RESULT,1);
NN=sum(CLASS==-1);
NP=sum(CLASS==1);

[RESULT_S,I]=sort(RESULT);
CLASS_S=CLASS(I);

TH=RESULT_S(NRSAM);
SAMNR=NRSAM;
TP=0;  TPa = [];
FP=0;  FPa = [];
TN=NN; TNa = [];
FN=NP; FNa = [];
TPR=0;
FPR=0;
AREA=0;
Q1B=0;
Q2B=0;
THRES_ROC=[TH];
TPR_ROC=[TPR];
FPR_ROC=[FPR];
SPEC_ROC=[TN/(FP+TN)];
ACC_ROC=[(TP+TN)/(NN+NP)];
PPV_ROC=[NaN];
NPV_ROC=[TN/(TN+FN)];


while ~isempty(TH)
   DELTA=CLASS_S(find(RESULT_S==TH));
   DFP=sum(DELTA==-1);
   DTP=sum(DELTA==1);
   TN=TN-DFP;
   AREA=AREA + DFP*TP + 0.5*DFP*DTP;
   Q2B=Q2B+DTP*((TN^2)+(TN*DFP)+((1/3)*(DFP^2)));
   Q1B=Q1B+DFP*((TP^2)+(TP*DTP)+((1/3)*(DTP^2)));
   FP=FP+DFP;
   TP=TP+DTP;
   FN=FN-DTP;
   TPR=TP/(TP+FN);
   FPR=FP/(FP+TN);
   
   SAMNR=max(find(RESULT_S<TH));
   TH=RESULT_S(SAMNR);
   
   TPR_ROC=[TPR_ROC ; TPR];
   FPR_ROC=[FPR_ROC ; FPR];
   THRES_ROC=[THRES_ROC ; TH];
   SPEC_ROC=[SPEC_ROC ; TN/(FP+TN)];
   ACC_ROC=[ACC_ROC ; (TP+TN)/(NN+NP)];
   if (TP+FP)==0
       PPV_ROC=[PPV_ROC ; NaN];
   else
       PPV_ROC=[PPV_ROC ; TP/(TP+FP)];
   end
   if (TN+FN)==0
       NPV_ROC=[NPV_ROC ; NaN];
   else
       NPV_ROC=[NPV_ROC ; TN/(TN+FN)];
   end

TPa = [TPa TP];
TNa = [TNa TN];
FPa = [FPa FP];
FNa = [FNa FN];
end

THRES_ROC=[THRES_ROC ; -1];

AREA=AREA/(NN*NP);
Q2=Q2B/((NN^2)*NP);
Q1=Q1B/(NN*(NP^2));

%Q1=AREA/(2-AREA);
%Q2=2*(AREA^2)/(1+AREA);
SE=sqrt((AREA*(1-AREA) + (NP-1)*(Q1-(AREA^2)) + (NN-1)*(Q2-(AREA^2)))/(NN*NP));



if fig(1)=='f',
  figure
  %fill([1 FPR_ROC 0],[0 TPR_ROC 0]','b');drawnow;
  plot(FPR_ROC,TPR_ROC,'b-','linewidth',2);
  title(['Receiver Operating Characteristic curve, area=' num2str(AREA) ...
       ', std = ',num2str(SE)]);
  xlabel('1 - Specificity');
  ylabel('Sensitivity');
end

FPR_ROC = FPR_ROC(2:end);
TPR_ROC = TPR_ROC(2:end);
