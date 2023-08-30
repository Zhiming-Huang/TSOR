function [Reward_summary,Vi_hat,PSi_p_hat] = TSOR(T,R,N,PSi,Obs,sigma_t,n_max,A,c_i,r_i)

T_index = zeros(T,1);
Loc_index = zeros(T,1);
Vi_hat = zeros(N,1);
Vi_hat(N) = -R;
Cost_summary = zeros(T,1);
Reward_summary = zeros(T,1);
PSi_p_hat = cell(N,1);
PSi_p_hat(N) = {1};
% Initialize the beta-distribution parameters
M_alpha = cell(N,1);
M_beta = cell(N,1);
for i = 1:N
    psi_cell = PSi(i);
    M_alpha(i) = {ones(1,numel(psi_cell{1,1}))};
    M_beta(i) = {ones(1,numel(psi_cell{1,1}))};
end

%%%At each time step
for i = 1:n_max
    % check if there are new packets?
    new_ps = find(sigma_t == i);
    % if not empty, than need to update value vector V_hat
    if ~isempty(new_ps)
        % sample the PSi_p_hat from beta distribution
        for j = 1:N-1
            PSi_p_hat(j) = {betarnd(cell2mat(M_alpha(j)),cell2mat(M_beta(j)))};
        end
        % Get the value of Vi_hat
        Vi_hat = calV(Vi_hat, A, PSi, PSi_p_hat, N, c_i, r_i);
        % Set the index
        T_index(new_ps) = 1;
        Loc_index(new_ps) = 1;
    end
    % Routing packets with the Vi_hat
    obs_n = Obs(:,:,i);
    current_packets = find(T_index == 1);
    for j = 1:numel(current_packets)
        packet_index = current_packets(j);
        current_location = Loc_index(packet_index);
        % Add one unit of cost
        Cost_summary(packet_index) = Cost_summary(packet_index) + c_i(current_location);
        % If exceed the Max R, drop it with 0 reward
        if Cost_summary(packet_index) > R
            T_index(packet_index) = 0;
        end
        % Observing the neighbors
        neighbours = find(obs_n(current_location,:) == 1);
        % update the M_alpha, M_beta
        alpha_v = cell2mat(M_alpha(current_location));
        beta_v = cell2mat(M_beta(current_location));
        psi_cell = PSi(current_location);
        psi_cell = psi_cell{1,1};
        neighbour_num = numel(alpha_v);
        if numel(neighbours) == 0
            alpha_v(1) = alpha_v(1) + 1;
            beta_v(2:neighbour_num) = beta_v(2:neighbour_num) + 1;
        else
            for k = 1:neighbour_num
                if isequal(cell2mat(psi_cell(k)), neighbours)
                    alpha_v(k) = alpha_v(k) + 1;
                else
                    beta_v(k) = beta_v(k) + 1;
                end
            end
        end
        M_alpha(current_location) = {alpha_v};
        M_beta(current_location) = {beta_v};
        % if the neighbours is not empty, then select the neighbour with
        % smallest value as the relay
        if ~isempty(neighbours)
            % if the destination is in the neighbours, then finish the
            % transmission, else send to the neighbours
            neighbours = [neighbours,current_location];
            [min_value, min_index] = min(Vi_hat(neighbours));
            if neighbours(min_index) == N
                T_index(packet_index) = 0;
                Loc_index(packet_index) = N;
                Reward_summary(packet_index) = R - Cost_summary(packet_index);
            else
                Loc_index(packet_index) = neighbours(min_index);
            end
        end
    end
end

disp('TSOR End');