function [msd]=getMSD(data,paramfile,delta_t,synID)

msd = struct;
run(paramfile)
global userParam
motion_data = data(synID).dat(:,1:2);
good_tp = size(motion_data,1);
 for ii=1:(good_tp-1) %lagtimes %
            dispacement2 = [];
            for jj=1:(good_tp-ii)% loop over all pairs separated by lagtimes
                % if use jj=1:ii:(good_tp-ii) then the
                % averaging will be over only independent pairs of points separated by lagtimes
                tp1 =jj;
                tp2 = jj+ii;
                dx=(motion_data(tp2,1)-motion_data(tp1,1))*userParam.pxtomicron;
                dy = (motion_data(tp2,2)-motion_data(tp1,2))*userParam.pxtomicron;
                % squared dispalcement traveled by cell celnter in time from tp1 to tp2
                d2 = ((power(dx,2)+power(dy,2)));%
                dt = ((tp2-tp1))*(delta_t/60);% time interval in hours
                TP(ii).times(jj,1:3) = [tp1 tp2 dt];
                dispacement2(jj,1:2) = [jj;  d2];% in micons
            end
            % disp(size(dispacement2(:,2),1));
            n_lags = (good_tp-ii);%size(dispacement2(:,2),1);
            %disp(((size(nonzeros(dispacement2(:,2)),1)*delta_t)/60))
            msd(synID).dat(ii,1) = sum((dispacement2(:,2)))/n_lags; %*delta_t)/60  [microns^2/hour] lag time = delta_t;
            % average over number of time steps at each lag time
            msd(synID).dat(ii,2) = ii; % how many delta_t intervals taken as the lag time
            msd(synID).trace_lengths = (good_tp);
            
            
 end  
end