function [Reward_summary,Vi_hat,PSi_p_hat,Value_learned_Node_2] = DORL(T,R,N,PSi,Obs,sigma_t,n_max,A,c_i,r_i,G)

T_index_explore = zeros(T,1);
T_index_exploit = zeros(T,1);
visitedlog = zeros(N,1);
Loc_index = zeros(T,1);
Dst_index = zeros(T,1);
n_mat = zeros(N, N);
m_mat = zeros(N, N);
o_i = zeros(N,1);
% l = N * ones(T,1);
Vi_hat = zeros(N,N);
Vi_hat(N,N) = -R;
Cost_summary = zeros(T,1);
Reward_summary = zeros(T,1);
Value_learned_Node_2 = zeros(T,N);
PSi_p_hat = cell(N,1);
number_explore_packets = 0;
number_total_packets = 0;
Indicator_psi = cell(N,1);
for i = 1:N
    psi_cell = PSi(i);
    Indicator_psi(i) = {zeros(1,numel(psi_cell{1,1}))};
    br = rand(1,numel(psi_cell{1,1}));
    br = br/sum(br);
    PSi_p_hat(i) = {br};
end

%%%At each time step
for i = 1:n_max
    % check if there are new packets?
    new_ps = find(sigma_t == i);
    % if not empty, than need to update value vector V_hat
    if ~isempty(new_ps)
        % Get the value of Vi_hat
%         Vi_hat = calV(Vi_hat, A, PSi, PSi_p_hat, N, c_i, r_i);
        % Set the index
%         T_index(new_ps) = 1;
        Loc_index(new_ps) = 1;
        for new_ps_index = 1:numel(new_ps)       
            number_total_packets = number_total_packets + 1;
            if number_explore_packets < (G * log(number_total_packets)) || number_total_packets == 1
                number_explore_packets = number_explore_packets + 1;
                T_index_explore(new_ps(new_ps_index)) = 1;
%                 if number_total_packets == 1
%                     Dst_index(new_ps(new_ps_index)) = N;
%                 else
                    % calculate the l(t - 1)
%                     l_minus_1 = min(n_i(number_total_packets-1, :));
%                     Dst_index(new_ps(new_ps_index)) = l_minus_1;
%                 end
            else
                T_index_exploit(new_ps(new_ps_index)) = 1;
            end
        end
    end
    obs_n = Obs(:,:,i);
    psi_update_flag = zeros(N,1);
    % Routing Exploration packets first
    current_explore_packets = find(T_index_explore == 1);
    for j = 1:numel(current_explore_packets)
        p_index = current_explore_packets(j);
        current_location = Loc_index(p_index);
        % Add one unit of COST
        Cost_summary(p_index) = Cost_summary(p_index) + c_i(current_location);
        % TTL expires
        if Cost_summary(p_index) > R
            T_index_explore(p_index) = 0;
            Cost_summary(p_index) = Cost_summary(p_index) - 1;
            Value_learned_Node_2(p_index,:) = Vi_hat(2,:);
            continue;
        end
        % Observing the neighbors
        neighbours = find(obs_n(current_location,:) == 1);
        allneighbours_index = find(A(current_location,:) == 1);
        allneighbours_index = allneighbours_index(allneighbours_index ~= N);
        % Add n_i
        n_mat(current_location, current_location) = n_mat(current_location, current_location) + 1;
        n_i_change_flag = 1;
        % The 1 case n_i increases by one, update m_i and o_i
        n_j = n_mat(current_location,allneighbours_index);
        [min_n_j_value, min_n_j_index] = min(n_j);
        original_m_i = m_mat(current_location, current_location);
        m_i_change_flag = 0;
        if original_m_i ~= min_n_j_value
            m_mat(current_location, current_location) = min_n_j_value;
            m_i_change_flag = 1;
        end
        o_i(current_location) = allneighbours_index(min_n_j_index);
        % Start to local broadcast, a decentralized way
        if ~isempty(neighbours)
            % each neighbours will receive the new n_j or m_j from
            % currentlocation
            for ix = 1:numel(neighbours)
                neighbor_index = neighbours(ix);
                % Use recursive function to implement
                [m_mat,n_mat,o_i] = update_neighbor(current_location,neighbor_index,n_i_change_flag,m_i_change_flag,m_mat,n_mat,A,obs_n,o_i);
            end
        end
        
        % Update the P_s_i_hat
        if psi_update_flag(current_location) == 0
            indicator_v = cell2mat(Indicator_psi(current_location));
            psi_cell = PSi(current_location);
            psi_cell = psi_cell{1,1};
            neighbour_num = numel(indicator_v);
            if numel(neighbours) == 0
                indicator_v(1) = indicator_v(1) + 1;
            else
                for k = 1:neighbour_num
                    if isequal(cell2mat(psi_cell(k)), neighbours)
                        indicator_v(k) = indicator_v(k) + 1;
                    end
                end
            end
            Indicator_psi(current_location) = {indicator_v};
            PSi_p_hat(current_location) = {indicator_v / n_mat(current_location, current_location)};
            psi_update_flag(current_location) = 1;
            % Update the vihat
            % Get the values from its neighbors
            if ~isempty(neighbours)
                for ix = 1:numel(neighbours)
                    neighbor_index = neighbours(ix);
                    Vi_hat(current_location, neighbor_index) = Vi_hat(neighbor_index, neighbor_index);
                end
            end
            allneighbours_index = find(A(current_location,:) == 1);
            Vi_current = Vi_hat(current_location, allneighbours_index);
            [Vi_current_value, Vi_current_index] = sort(Vi_current);
            % Rank neighbours with Vi
            allneighbours_index = allneighbours_index(Vi_current_index);
            for rank_s = 1:numel(allneighbours_index) + 1
                Set_A_flag = zeros(N,1);
                Set_X_flag = zeros(N,1);
                Set_A_flag(allneighbours_index(1:rank_s-1)) = 1;
                Set_X_flag([allneighbours_index,current_location]) = 1;
                Set_X_flag(allneighbours_index(1:rank_s-1)) = 0;
                Vi = Vi_hat(current_location, :);
                Vi = calV_de(Vi, A, PSi, PSi_p_hat, c_i, r_i, Set_A_flag, Set_X_flag);
                [Vi_value, Vi_index] = sort(Vi);
                real_rank = find(Vi_index == current_location);
                if real_rank == rank_s
                    Vi_hat(current_location, current_location) = Vi(current_location);
                    break;
                end
            end
        end
        % If the neighbours are not empty
        if ~isempty(neighbours)
            % not reach l(t-1) node yet
            for ix = 1:numel(neighbours)
                neighbor_index = neighbours(ix);
                if neighbor_index == N
                    T_index_explore(p_index) = 2;
                    Loc_index(p_index) = N;
                    Reward_summary(p_index) = R - Cost_summary(p_index);
                    Value_learned_Node_2(p_index,:) = Vi_hat(2,:);
                    continue;
                elseif neighbor_index == o_i(current_location)
%                     T_index_explore(p_index) = o_i(current_location);
                    Loc_index(p_index) = o_i(current_location);
                    visitedlog(current_location) = visitedlog(current_location)+1;
                end
            end
            % reached l(t-1) node
%             else
%                 neighbours = [neighbours,current_location];
%                 [min_value, min_index] = min(Vi_hat(current_location, neighbours));
%                 if neighbours(min_index) == N
%                     T_index_exploit(packet_index) = 2;
%                     Loc_index(packet_index) = N;
%                     Reward_summary(packet_index) = R - Cost_summary(packet_index);
%                 else
%                     Loc_index(packet_index) = neighbours(min_index);
%                 end
%             end
        end
    end
    
    % Routing Exploitation packets
    current_exploit_packets = find(T_index_exploit == 1);
    for j = 1:numel(current_exploit_packets)
        packet_index = current_exploit_packets(j);
        current_location = Loc_index(packet_index);
        % Add one unit of COST
        Cost_summary(packet_index) = Cost_summary(packet_index) + c_i(current_location);
        % TTL expires
        if Cost_summary(packet_index) > R
            T_index_exploit(packet_index) = 0;
            Cost_summary(packet_index) = Cost_summary(packet_index) - 1;
            Value_learned_Node_2(packet_index,:) = Vi_hat(2,:);
            continue;
        end
        neighbours = find(obs_n(current_location,:) == 1);
%         allneighbours_index = find(A(current_location,:) == 1);
%         % Add n_i
%         n_mat(current_location, current_location) = n_mat(current_location, current_location) + 1;
%         n_i_change_flag = 1;
%         % The 1 case n_i increases by one, update m_i and o_i
%         n_j = n_mat(current_location,allneighbours_index);
%         [min_n_j_value, min_n_j_index] = min(n_j);
%         original_m_i = m_mat(current_location, current_location);
%         m_i_change_flag = 0;
%         if original_m_i ~= min_n_j_value
%             m_mat(current_location, current_location) = min_n_j_value;
%             m_i_change_flag = 1;
%         end
%         o_i(current_location) = allneighbours_index(min_n_j_index);
%         % Start to local broadcast, a decentralized way
%         if ~isempty(neighbours)
%             % each neighbours will receive the new n_j or m_j from
%             % currentlocation
%             for ix = 1:numel(neighbours)
%                 neighbor_index = neighbours(ix);
%                 % Use recursive function to implement
%                 [m_mat,n_mat,o_i] = update_neighbor(current_location,neighbor_index,n_i_change_flag,m_i_change_flag,m_mat,n_mat,A,obs_n,o_i);
%             end
%         end
        
        % if the destination is in the neighbours, then finish the
        % transmission, else send to the neighbours
        neighbours = [neighbours,current_location];
        [min_value, min_index] = min(Vi_hat(current_location, neighbours));
        if neighbours(min_index) == N
            T_index_exploit(packet_index) = 2;
            Loc_index(packet_index) = N;
            Reward_summary(packet_index) = R - Cost_summary(packet_index);
            Value_learned_Node_2(packet_index,:) = Vi_hat(2,:);
        else
            Loc_index(packet_index) = neighbours(min_index);
        end
    end
end

% disp('DORL End');