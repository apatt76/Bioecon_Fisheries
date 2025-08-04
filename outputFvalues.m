function status = outputFvalues(t, y, flag, params)
    persistent F_list
    persistent J_list

    if isempty(flag)
        F = computeF(y, params);
        F_list(end+1,:,:) = F;

        J = computeJacobian(@differential, y, t, params);
        J_list(end+1,:,:) = J;

    elseif strcmp(flag, 'init')
        F_list = [];
        J_list = [];

    elseif strcmp(flag, 'done')
        assignin('base', 'F_series', F_list);
        assignin('base', 'J_series', J_list);
    end

    status = 0;
end
