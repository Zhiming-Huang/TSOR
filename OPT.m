function Reward_summary = OPT(Vi_star,T,R,N,Obs,sigma_t,n_max,c_i)

T_index = zeros(T,1);
Loc_index = zeros(T,1);
Cost_summary = zeros(T,1);
Reward_summary = zeros(T,1);

%%%At each time step
for i = 1:n_max
    % check if there are new packets?
    new_ps = find(sigma_t == i);
    % if not empty, than need to update value vector V_hat
    if ~isempty(new_ps)
        % Set the index
        T_index(new_ps) = 1;
        Loc_index(new_ps) = 1;
    end
    % Routing packets with the Vi_star
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
        % if the neighbours is not empty, then select the neighbour with
        % smallest value as the relay
        if ~isempty(neighbours)
            % if the destination is in the neighbours, then finish the
            % transmission, else send to the neighbours
            neighbours = [neighbours, current_location];
            [min_value, min_index] = min(Vi_star(neighbours));
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

disp('OPT end');