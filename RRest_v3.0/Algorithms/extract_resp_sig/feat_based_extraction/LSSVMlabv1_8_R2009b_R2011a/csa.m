function [pfinal,efinal] = csa(pn,herrfunc,varargin)

%
% Internal function based on 
% Xavier-de-Souza S, Suykens JA, Vandewalle J, Bolle D., Coupled Simulated
% Annealing, IEEE Trans Syst Man Cybern B Cybern. 2010 Apr;40(2):320-35.
%
% Copyright (c) 2011,  KULeuven-ESAT-SCD, License & help @% http://www.esat.kuleuven.be/sista/lssvmlab


%switch length(varargin)
OPT.T0 = 1; OPT.Tac0 = 1; OPT.FEmax = 40; OPT.FTsteps = 20; OPT.etol = 1e-45; OPT.print = 0;
T0 = OPT.T0;  % initial temperature
Tac0 = OPT.Tac0;  % initial temperature
FEmax = OPT.FEmax; % max # of function evaluations
FTsteps = OPT.FTsteps; % # of steps at fix temperature
etol = OPT.etol; % energy tolerance

clear OPT
% initializes M
pdim = size(pn,1);
pnum = size(pn,2);

NT = ceil(FEmax/FTsteps/pnum); % # max. number of cooling cycles
NI = FTsteps; % #steps per temperature

rand('twister',sum(pn(1)*clock))
randn('state',sum(pn(2)*clock))

e0 = feval(herrfunc,pn,varargin{:});
%if any(e0<0), etol = -etol*etol^-2;end

p0 = pn;
[be0,ind] = min(e0);
bp0 = pn(:,ind);

pblty = zeros(1,pnum);
sgnCR = -1;
CR = 0.1;%0.05;
pvar_est = 0.995;

Tac = Tac0;

for k = 1:NT,
    C = progress('init',['Determine initial tuning parameters for simplex...',': # cooling cycle(s) ', num2str(k)]);
    pbltvar = var(pblty,1);
    sgnCR_ant = sgnCR;
    sgnCR = 2*((pbltvar > (pvar_est*(pnum-1)/(pnum^2)))-0.5);
    Tac = Tac + sgnCR*CR*Tac;

    % T schedules
    % T = T0/log(k+1);
    T = T0/k;
    %    T = T0*exp(-k);
    %    T = T0*exp(-exp(k));

    % Tac schedules (not needed when variance control is applied)
    %    Tac = Tac0/log(h+1);
    %    Tac = Tac0/h;
    %    Tac = Tac0*exp(-k);
    %    Tac = Tac0*exp(-exp(k));

    %    TV(k) = Tac;
    for l = 1:NI,
        % choose new coordinates and compute
        % the function to minimize
        r = tan(pi*(rand(pdim,pnum)-0.5));
        % ****************** Wrapping ***************
        %pn = 2*mod((p0 + r * T + 1)/2,1)-1;
        %**************** Non wrapping *************
         pn = p0 + r * T ;%* (diag(e0)./sum(e0));
         indd = find(abs(pn)>10);
         while(numel(indd)),
             r(indd) = tan(pi*(rand(size(indd))-0.5));
             pn = p0 + r * T ;%* (diag(e0)./sum(e0));
             indd = find(abs(pn)>15);
         end

        en = feval(herrfunc,pn,varargin{:});
        Esum = sum(exp((e0-max(e0))./Tac));
        for i=1:pnum
            pblty(i) = min(1,exp((e0(i)-max(e0))/Tac)/Esum);
            if (en(i)-e0(i)) < 0,
                % accept
                p0(:,i) = pn(:,i); e0(i) = en(i);
                if e0(i) < be0, 
                    be0 = e0(i); 
                    bp0 = p0(:,i); 
%                     if OPT.print == 1,
%                         fprintf('v%e %d\n',min(e0),k);
%                     end
                end
            else
                r = rand;
                if (pblty(i)) >= r
                    % accept
                    p0(:,i) = pn(:,i); e0(i) = en(i);
                end
            end
        end
         C = progress(C,l/NI);
        if any(e0<etol), break, end
    end
    if any(e0<etol), break, end
end
% kfinal = ((k-1)*NI + l)*pnum;
efinal = be0;
pfinal = bp0;
