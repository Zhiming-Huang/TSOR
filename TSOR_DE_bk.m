function [Reward_summary,Vi_hat,PSi_p_hat] = TSOR_DE(T,R,N,PSi,Obs,sigma_t,n_max,A,A_p,c_i,r_i)

T_index = zeros(T,1);
Loc_index = zeros(T,1);
UP_index = zeros(T+5,1);
UP_index(1) = 1;
up_num = 1;
UP_Loc = zeros(T+5,N);
UP_Loc(1,N) = 1;
UP_TTL = R * ones(T+5,1);
Vi_hat = zeros(N,N);
Vi_hat(N,N) = -R;
Cost_summary = zeros(T,1);
Reward_summary = zeros(T,1);
PSi_p_hat = cell(N,1);
% PSi_p_hat(N) = {1};
% Initialize the beta-distribution parameters
M_alpha = cell(N,1);
M_beta = cell(N,1);
for i = 1:N
    psi_cell = PSi(i);
    M_alpha(i) = {ones(1,numel(psi_cell{1,1}))};
    M_beta(i) = {ones(1,numel(psi_cell{1,1}))};
%     Vi_hat(i) = {ones(1,numel(psi_cell{1,1}))};
end

% sample the PSi_p_hat from beta distribution
for j = 1:N
    PSi_p_hat(j) = {betarnd(cell2mat(M_alpha(j)),cell2mat(M_beta(j)))};
end
%% At each time step
for i = 1:n_max
    % check if there are new packets?
    new_ps = find(sigma_t == i);
    % if not empty, than need to update value vector V_hat
    if ~isempty(new_ps)
        % sample the PSi_p_hat from beta distribution
%         for j = 1:N
%             PSi_p_hat(j) = {betarnd(cell2mat(M_alpha(j)),cell2mat(M_beta(j)))};
%         end
        % Set the index
        T_index(new_ps) = 1;
        Loc_index(new_ps) = 1;
    end
    % Current link status
    obs_n = Obs(:,:,i);
    psi_update_flag = zeros(N,1);
    %% Routing regular packets with the Vi_hat
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
        % update the M_alpha, M_beta, if not updated in this slot
        if psi_update_flag(current_location) == 0
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
            psi_update_flag(current_location) = 1;
            PSi_p_hat(current_location) = {betarnd(cell2mat(M_alpha(current_location)),cell2mat(M_beta(current_location)))};
        end
        % if the neighbours is not empty, then select the neighbour with
        % smallest value as the relay
        if ~isempty(neighbours)
            % if the destination is in the neighbours, then finish the
            % transmission, else send to the neighbours
            neighbours = [neighbours,current_location];
            [min_value, min_index] = min(Vi_hat(current_location, neighbours));
            if neighbours(min_index) == N
                T_index(packet_index) = 0;
                Loc_index(packet_index) = N;
                Reward_summary(packet_index) = R - Cost_summary(packet_index);
                up_num = up_num + 1;
                UP_index(up_num) = 1;
                UP_Loc(up_num, N) = 1;
            else
                Loc_index(packet_index) = neighbours(min_index);
            end
        end
    end
    %% ACK with value updated
    current_ups = find(UP_index == 1);
    for j = 1:numel(current_ups)
        up_index = current_ups(j);
        UP_TTL(up_index) = UP_TTL(up_index) - 1;
        % TTL timeout in next slot
        if UP_TTL(up_index) == 0
            UP_index(up_index) = 0;
        end
        % BROADCAST PROCESSINGGGGG
        % Before processing
        % Change the status of 2 to 1 in last round to keep them
        % broadcasting
        % change_index = find(UP_Loc(up_index,:) == 2);
        temp_locs = UP_Loc(up_index,:);
        if all(temp_locs == 3)
            UP_index(up_index) = 0;
            continue;
        end
        temp_locs(find(UP_Loc(up_index,:) == 2)) = 1;
        UP_Loc(up_index,:) = temp_locs;
        while 1
            locs = UP_Loc(up_index,:);
            broadcast_locs = find(locs == 1);
            if isempty(broadcast_locs)
                break;
            else
                for k = 1:numel(broadcast_locs)
                    current_location = broadcast_locs(k);
                    % Observing the neighbors
                    neighbours = find(obs_n(current_location,:) == 1);
                    % all neighbours  
                    allneighbours = find(A(current_location,:) == 1);
                    % update the M_alpha, M_beta, if not updated in this slot
                    if psi_update_flag(current_location) == 0
                        alpha_v = cell2mat(M_alpha(current_location));
                        beta_v = cell2mat(M_beta(current_location));
                        psi_cell = PSi(current_location);
                        psi_cell = psi_cell{1,1};
                        neighbour_num = numel(alpha_v);
                        if numel(neighbours) == 0
                            alpha_v(1) = alpha_v(1) + 1;
                            beta_v(2:neighbour_num) = beta_v(2:neighbour_num) + 1;
                        else
                            for l = 1:neighbour_num
                                if isequal(cell2mat(psi_cell(l)), neighbours)
                                    alpha_v(l) = alpha_v(l) + 1;
                                else
                                    beta_v(l) = beta_v(l) + 1;
                                end
                            end
                        end
                        M_alpha(current_location) = {alpha_v};
                        M_beta(current_location) = {beta_v};
                        psi_update_flag(current_location) = 1;
                        PSi_p_hat(current_location) = {betarnd(cell2mat(M_alpha(current_location)),cell2mat(M_beta(current_location)))};
                    end
                    if ~isempty(neighbours)
                        % For each neighbors
                        for ix = 1:numel(neighbours)
                            neighbor_index = neighbours(ix);
                            % if it already receives the broadcast packet,
                            % then only update the new calculated value,
                            % else add into the broadcast nodes and
                            % recalculate its updated value
                            Vi_hat(neighbor_index, current_location) = Vi_hat(current_location, current_location);
                            locs = UP_Loc(up_index,:);
                            if locs(neighbor_index) == 0
                                UP_Loc(up_index,neighbor_index) = 1;
                                % Calculate its new value
                                % neighbour_neighbours = find(obs_n(neighbor_index,:) == 1);
%                                 calV_de(Vi, A, PSi, PSi_p, c_i, r_i, Set_A_flag, Set_X_flag);
                                % Get the number of neighbors
                                neighbour_neighbours = find(A(neighbor_index,:) == 1);

                                %%V3
                                Vi_current = Vi_hat(neighbor_index, neighbour_neighbours);
                                [Vi_current_value, Vi_current_index] = sort(Vi_current);
                                % Rank neighbours with Vi
                                neighbour_neighbours = neighbour_neighbours(Vi_current_index);
                                for rank_s = 1:numel(neighbour_neighbours) + 1
                                    Set_A_flag = zeros(N,1);
                                    Set_X_flag = zeros(N,1);
                                    Set_A_flag(neighbour_neighbours(1:rank_s-1)) = 1;
                                    Set_X_flag([neighbour_neighbours,neighbor_index]) = 1;
                                    Set_X_flag(neighbour_neighbours(1:rank_s-1)) = 0;
                                    Vi = Vi_hat(neighbor_index, :);
                                    Vi = calV_de(Vi, A, PSi, PSi_p_hat, c_i, r_i, Set_A_flag, Set_X_flag);
                                    [Vi_value, Vi_index] = sort(Vi);
                                    real_rank = find(Vi_index == neighbor_index);
                                    if real_rank == rank_s
                                        Vi_hat(neighbor_index, neighbor_index) = Vi(neighbor_index);
                                        break;
                                    end
                                end
%% V2
%                                 Vi_current = Vi_hat(neighbor_index, neighbour_neighbours);
%                                 [Vi_current_value, Vi_current_index] = sort(Vi_current);
%                                 % Rank neighbours with Vi
%                                 neighbour_neighbours = neighbour_neighbours(Vi_current_index);
%                                 for rank_s = 1:numel(neighbour_neighbours) + 1
%                                     Set_A_flag = zeros(N,1);
%                                     Set_X_flag = zeros(N,1);
%                                     Set_A_flag(neighbour_neighbours(1:rank_s-1)) = 1;
%                                     Set_X_flag([neighbour_neighbours,neighbor_index]) = 1;
%                                     Set_X_flag(neighbour_neighbours(1:rank_s-1)) = 0;
%                                     Vi = Vi_hat(neighbor_index, :);
%                                     Vi = calV_de(Vi, A, PSi, PSi_p_hat, c_i, r_i, Set_A_flag, Set_X_flag);
%                                     [Vi_value, Vi_index] = sort(Vi);
%                                     real_rank = find(Vi_index == neighbor_index);
%                                     if real_rank == rank_s
%                                         Vi_hat(neighbor_index, neighbor_index) = Vi(neighbor_index);
%                                         break;
%                                     end
%                                 end
%% V1
%                                 Vi_current = Vi_hat(neighbor_index, neighbour_neighbours);
%                                 [Vi_current_value, Vi_current_index] = sort(Vi_current);
%                                 % Rank neighbours with Vi
%                                 neighbour_neighbours = neighbour_neighbours(Vi_current_index);
%                                 for rank_s = 1:numel(neighbour_neighbours) + 1
%                                     Set_A_flag = zeros(N,1);
%                                     Set_X_flag = zeros(N,1);
%                                     Set_A_flag(neighbour_neighbours(1:rank_s-1)) = 1;
%                                     Set_A_flag(neighbor_index) = 1;
%                                     Set_X_flag(neighbour_neighbours) = 1;
%                                     Set_X_flag(neighbour_neighbours(1:rank_s-1)) = 0;
%                                     Vi = Vi_hat(neighbor_index, :);
%                                     Vi = calV_de(Vi, A, PSi, PSi_p_hat, c_i, r_i, Set_A_flag, Set_X_flag);
%                                     [Vi_value, Vi_index] = sort(Vi);
%                                     real_rank = find(Vi_index == neighbor_index);
%                                     if real_rank == rank_s
%                                         Vi_hat(neighbor_index, neighbor_index) = Vi(neighbor_index);
%                                         break;
%                                     end
%                                 end
                            end
                        end
                    end
                    % temporary stop the broadcast at this round
                    UP_Loc(up_index,current_location) = 2;
                    % If all neighbours are already broadcasted, then stop
                    % broadcasting the packet in current location
                    temp_locs = UP_Loc(up_index,:);
                    if all(temp_locs(allneighbours) > 0)
                        UP_Loc(up_index,current_location) = 3;
                    end
                end
            end
        end
    end
end

if sum(Reward_summary) < 6000
end
disp('TSORDE end');