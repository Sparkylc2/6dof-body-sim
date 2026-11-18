function after_run(model, varargin)
%AFTER_RUN open bsim.rocket.replay once the model has completely finished running
%   bsim.rocket.after_run(bdroot, xlog, terms, Prk, 'U', ulog)      (in the model's StopFcn)
%
%   the StopFcn fires as soon as the solver's done, while simulink is still wrapping up (logging,
%   the run button still says stop), so this waits until the model reports it's fully stopped and
%   only then opens and plays the replay. everything after model goes straight to the replay
tm = timer('ExecutionMode', 'fixedSpacing', 'Period', 0.2, 'TasksToExecute', 300, ...
    'TimerFcn', @check, 'StopFcn', @(t, ~) delete(t));
start(tm);

    function check(t, ~)
        try
            status = get_param(model, 'SimulationStatus');
        catch
            stop(t); % the model's been closed
            return
        end
        if strcmp(status, 'stopped')
            stop(t);
            bsim.rocket.replay(varargin{:});
        end
    end
end
