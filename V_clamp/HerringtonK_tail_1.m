function [v_half_act_IK, n_act_IK] = HerringtonK_tail_1(parameters, parameter_names, states, state_names, rhs_model)
    % Herrington protocols
%     parameters(find(strcmp(parameter_names, 'g_KATP0'))) = 0.01;
%     parameters(find(strcmp(parameter_names, 'g_KATP_hat'))) = 0.05;
%     parameters(find(strcmp(parameter_names, 'g_SK'))) = 0.1;
%     parameters(find(strcmp(parameter_names, 'g_BK'))) = 0.02;
%     parameters(find(strcmp(parameter_names, 'g_Kv'))) = 1.0;
%     parameters(find(strcmp(parameter_names, 'g_HERG'))) = 0;
    parameters(find(strcmp(parameter_names, 'V_K'))) = -92;
    
    v_rest = -80;
    v_pulse = (-70:10:40);
    v_tail = -40;
    rest_length = 5000;
    pulse_length = 100;
    prepulse_length = 0;
    end_length = 0;
    tail_length = 20;
    num_pulses = length(v_pulse);
    cycle_length = rest_length+pulse_length+tail_length;
    current = 'K_tail_1';
    plotting = 'true';
    options = odeset('MaxStep',1);

    v_clamp_amps(1:3:num_pulses*3) = v_rest; v_clamp_amps(2:3:num_pulses*3) = v_pulse; v_clamp_amps(3:3:num_pulses*3) = v_tail;
    v_clamp_times(1:3:3*num_pulses) = rest_length:cycle_length:num_pulses*cycle_length; v_clamp_times(2:3:num_pulses*3) = rest_length+pulse_length:cycle_length:num_pulses*cycle_length; v_clamp_times(3:3:num_pulses*3)=cycle_length:cycle_length:num_pulses*cycle_length; 
    R_clamp = 0.001; % (GOhm)
  
    init_tspan = [0,10000];
    [T_init, Y_init] = ode15s(rhs_model,init_tspan, states, options, parameters, v_rest, init_tspan(end), R_clamp);
    
    tspan = [0,v_clamp_times(end)];
    [T, Y] = ode15s(rhs_model,tspan, Y_init(end,:), options, parameters, v_clamp_amps, v_clamp_times, R_clamp);

    try
        outs = analyze_v_clamp(Y,T,parameters,state_names, cycle_length, prepulse_length, end_length, rhs_model, v_clamp_amps, v_clamp_times, R_clamp, current, plotting);
        v_half_act_IK = outs{1};
        n_act_IK = outs{2};
    catch ME
        reason = {ME.identifier};
        disp(reason)
        v_half_act_IK = NaN;
        n_act_IK = NaN;
    end
end