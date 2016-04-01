function rr = TFu(data, up)
%TFu fuses RR estimates using temporal fusion.

beta = 0.8;

rr.t = data.t;
rr.v = nan(length(data.v),1);

% if there was actually an RR estimated:
if ~(sum(isnan(data.v)) == length(data.v))
    
    rr.v(1) = data.v(min(find(~isnan(data.v))));
    for win_no = 2 : length(rr.t)
        if ~isnan(data.v(win_no))
            rr.v(win_no) = ( (1-beta)*data.v(win_no) ) + ( beta*rr.v(win_no-1) );
        else
            rr.v(win_no) = rr.v(win_no-1);
        end
    end

end


end