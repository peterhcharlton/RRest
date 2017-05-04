function [xopt,fopt]=simann(func, x, LB, UB, sa_t, sa_rt, sa_nt, sa_ns,rseed) 

% Simulated Annealing programmed for minimization problem
% INPUTS
% func, string variable containing name of function file to be optimized 
% x, starting values
% LB, lower bound on optimization parameters
% UB, upper bound on optimization parameters
% sa_t, initial temperature
% sa_rt, temperature reduction factor, 0 < sa_rt < 1, try .85
% sa_nt, number of times through ns loop before temperature reduction (recommended value: 5)
% sa_ns, number of times through function before stepsize adjustment (recommended value: 20)
%
% OUTPUTS
% xopt, the optimal solution
%
% 


LB=LB(:)';                                   
UB=UB(:)';

rand('state',rseed);                      %sets seed for random number generator
sa_neps=4;                                %number of times eps
                                          %tolerance is achieved before termination
sa_eps=eps;                              %convergence criteria
sa_maxeval=60;%12000000;                      %maximum number of function evaluations

sa_nargs=length(LB);                      %number of parameters
sa_nobds=0;
sa_nacc=0;                                %number of acceptions
sa_nevals=0;                              %number of evaluations
sa_opteval=0;                             %optimum number of
                                          %function evaluations

fstar=Inf*ones(sa_neps,1);

%x=LB+(UB-LB).*rand(1, sa_nargs);         %starting values for model parameters
f=feval(func,x);                          %function evaluation with parameters x
%disp('initial loss function value:');disp(f);
sa_nevals=sa_nevals+1;
xopt=x;
fopt=f;
xtot=x;
fstar(1)=f;

VM=(UB-LB);%/2;                      %maximum step size

%LOOP
while 1 
  
  nup=0;                                   %number of uphill movements
  nrej=0;                                  %number of rejections
  nnew=0;
  ndown=0;                                 %number of downhill movements
  lnobds=0;
  nacp=zeros(sa_nargs,1);
  C = progress('init','Determine initial hyperparameters for simplex...');
  for m=1:sa_nt
    for j=1:sa_ns
      for h=1:sa_nargs
        if sa_nevals>=sa_maxeval
          %disp('too many function evaluations')
          return
        end
        C = progress(C,sa_nevals/sa_maxeval);
        %workbar(sa_nevals/sa_maxeval,'Determine initial hyperparameters to build grid...','Progress') 
        % generate xp, trial value of x
        xp=x;
        xp(h)=x(h)+VM(h)*(2*rand(1,1)-1.0);               %calculate new value for x (xp)
        if (xp(h)<LB(h)) | (xp(h)>UB(h))
          xp(h)=LB(h)+(UB(h)-LB(h))*rand(1,1);
          lnobds=lnobds+1;
          sa_nobds=sa_nobds+1;
        end       
        % evaluate at xp and return as fp
        %disp ('current parameter vector:');disp(xp);
        fp=feval(func,xp);                                 %function evaluation with parameters xp
        %disp ('function value');disp(fp);
        sa_nevals=sa_nevals+1;

        % we minimize! accept if the function value decreases
        if fp<=f
          x=xp;
          f=fp;
          sa_nacc=sa_nacc+1;
          nacp(h)=nacp(h)+1;
          nup=nup+1;
          % if smaller than any previous point, record as new optimum
          if fp<fopt
            xopt=xp;
            fopt=fp;
            sa_opteval=sa_nevals;
            nnew=nnew+1;
          end
        else % function value increases
          p=exp((f-fp)/sa_t);                              %random number
          pp=rand(1,1);
          if pp<p
            x=xp;
            f=fp;
            sa_nacc=sa_nacc+1;
            nacp(h)=nacp(h)+1;
            ndown=ndown+1;
          else
            nrej=nrej+1;
          end
        end
      end
    end
    
    % adjust maximal step size vm
    c=ones(sa_nargs,1)*2;%??
    for i=1:sa_nargs
      ratio=nacp(i)/sa_ns;
      if ratio>0.6
        VM(i)=VM(i) * (1+c(i)*(ratio-0.6)/0.4);
      elseif ratio <0.4
        VM(i)=VM(i)/(1+c(i)*((0.4-ratio)/0.4));
      end
      if VM(i)>(UB(i)-LB(i))
        VM(i)=UB(i)-LB(i);
      end
    end

    % provide statistics about current state of optimization
    
%     disp('No. of evaluations');disp(sa_nevals);disp('  current temperature');disp(sa_t);
%     disp('current optimum function value');disp(fopt);
%     disp('No. of downhill steps');disp(nup);  % note misnomer in variable declaration!
%     disp('No. of accepted uphill steps');disp(ndown); % we minimize, thus downhill is always accepted!
%     disp('No. of rejections');disp(nrej);
%     disp('current parameter values');disp(xp);
%     disp('current optimum vector');disp(xopt);
%     disp('current step size');disp(VM);
    %disp('Variables used:');whos;

  for i=1:sa_nargs
     nacp(i) = 0;
  end
  end
  
  
  % check termination criteria
  fstar(1)=f;
  quit = ((fstar(1)-fopt) <= sa_eps);
  if any(abs(fstar-f)>sa_eps)
    quit=0;
  end
  
  if quit
    disp(['simulated annealing achieved termination after ', num2str(sa_nevals),' evals']);
    return
  end
  
  % reduce temperature  
  sa_t=sa_t*sa_rt;
  fstar(2:4)=fstar(1:3);
  % continue from current optimum
  x=xopt;
  f=fopt;
end %while



