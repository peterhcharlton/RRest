function cont=progress(status,rate)
%PROGRESS   Text progress bar
%   Similar to waitbar but without the figure display.
%
%   Start:
%      C = PROGRESS('init',TITLE), where the default TITLE is 'please wait'
%
%   On progress:
%      C = PROGRESS(C,RATE);
%
%   Examples:
%      n=20;
%      for i=1:n
%        if i==1, c=progress('init');
%        else     c=progress(c,i/n); end
%
%        % computing something ...
%        pause(.1)
%      end
%
%      % inside a script you may use:
%       c=progress('init','wait for ... whatever');
%       for i=1:n
%         c=progress(c,i/n);
%         ...
%       end
%

lmax=50;
if isequal(status,'init')
  if nargin > 1
    title = rate;
    fprintf(1,'\n   %s\n',title);
    str = repmat(' ',1,lmax-4);
  else
    str = '                please wait                   ';
    fprintf(1,'\n');
  end
  fprintf(1,'  |-%s-|\n',str);
  fprintf(1,'%s','  ');
  cont=0;
else
  cont = status;
  n=ceil(rate*lmax);
  dif = n - cont;
  if dif > 0 & n <=lmax
    N=dif;
  else
    N=0;
  end
  cont=cont+N;
  str = repmat('*',1,N);
  fprintf(1,'%s',str);
  if rate==1
   fprintf(1,'%s\n',' done');
  end
end

