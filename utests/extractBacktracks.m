function [logic_backtracks, varargout] = extrBacktracks(time,pwcs)
%EXTRORIGBACKTRACKS returns a logic array were intervals of the transcription trace
%considered as backtracking are true and false otherwise
%   Detailed explanation goes here

    %transpose arrays if necessary
    if(size(pwcs,1)>size(pwcs,2))
        pwcs = pwcs';
    end
    if (pwcs(1) ~= pwcs(2)) % trace shoudn't start with a step
        pwcs(1) = pwcs(2);  % since there's no corresponding dwell time
    end
    if(size(time,1)>size(time,2))
        time = time';
    end

    %separate forward steps of elongation from backtracking steps
    x_forward = pwcs;
    for i=2:length(pwcs)
            if(pwcs(i)-x_forward(i-1)<=0) %(pwcs(i)-x_forward(i-1)<=0)
                x_forward(i) = x_forward(i-1);
            end
    end
    logic_backtracks = (pwcs-x_forward)<0;
    
    posIndex = 1:length(pwcs);
    %step_pts = [pwcs(diff([pwcs(1) pwcs])~=0)];
    %step_times = [time(diff(pwcs)~=0)];
    %array_pts = [1 posIndex(diff(pwcs)~=0)];
    
    stepPos_logic = diff([pwcs(1) pwcs])~=0;
    stepPos_all = posIndex(stepPos_logic~=0);
    stepPos_backtr = posIndex(logic_backtracks~=0);
    
    for i=1:length(stepPos_backtr)
        
        ind = find(stepPos_all==stepPos_backtr(i));
        if(ind>1)
            logic_backtracks(stepPos_all(ind-1):stepPos_all(ind)) = 1;
        end
    end
    
    
    s=0;
    if(logic_backtracks(1))
        s=s+1;
        backtr_pos{s} = pwcs(1);
        backtr_ti{s} = time(1);
    end
    for i=2:length(logic_backtracks)
        if(logic_backtracks(i))
            if(~logic_backtracks(i-1))
                s = s+1;
                backtr_pos{s} = [];
                backtr_ti{s} = [];
            end
            backtr_pos{s} = [backtr_pos{s} pwcs(i)];
            backtr_ti{s} = [backtr_ti{s} time(i)];
        end
    end
    
    if nargout>1
        switch nargout
            
            case 2
               varargout{1} = backtr_ti; 
            case 3
                varargout{1} = backtr_ti;
                varargout{2} = backtr_pos;
        end
    end
end



