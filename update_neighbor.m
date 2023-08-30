function [m_mat,n_mat,o_i] = update_neighbor(current_location,neighbor_index,n_i_change_flag,m_i_change_flag,m_mat,n_mat,A,obs_n,o_i)

neighbor_neighbors = find(obs_n(neighbor_index,:) == 1);
neighbor_allneighbours_index = find(A(neighbor_index,:) == 1);
neighbor_allneighbours_index = neighbor_allneighbours_index(neighbor_allneighbours_index ~= numel(o_i));
% Case 2 or Case 3
if (n_i_change_flag && (n_mat(neighbor_index, current_location) < n_mat(current_location, current_location))) || (m_i_change_flag && (m_mat(neighbor_index, current_location) < m_mat(current_location, current_location)))
    % update the n_i in n_j
    if n_i_change_flag
        n_mat(neighbor_index, current_location) = n_mat(current_location, current_location);
    end
    % update the m_i in m_j
    if m_i_change_flag
        m_mat(neighbor_index, current_location) = m_mat(current_location, current_location);
    end
    n_j = n_mat(neighbor_index,neighbor_allneighbours_index);
    [min_n_j_value, min_n_j_index] = min(n_j);
    original_m_i = m_mat(neighbor_index, neighbor_index);
    neighbor_m_i_change_flag = 0;
    if original_m_i ~= min_n_j_value
        m_mat(neighbor_index, neighbor_index) = min_n_j_value;
        neighbor_m_i_change_flag = 1;
    end
    o_i(neighbor_index) = neighbor_allneighbours_index(min_n_j_index);
    if ~isempty(neighbor_neighbors)
        for ix = 1:numel(neighbor_neighbors)
            neighbor_neighbor_index = neighbor_neighbors(ix);
            % Use recursive function to implement
            [m_mat,n_mat,o_i] = update_neighbor(neighbor_index,neighbor_neighbor_index,0,neighbor_m_i_change_flag,m_mat,n_mat,A,obs_n,o_i);
        end
    end
elseif m_i_change_flag && (m_mat(neighbor_index, current_location) > m_mat(current_location, current_location)) && (m_mat(neighbor_index, neighbor_index) > m_mat(current_location, current_location))
    m_mat(neighbor_index, current_location) = m_mat(current_location, current_location);
    m_mat(neighbor_index, neighbor_index) = m_mat(current_location, current_location);
    o_i(neighbor_index) = current_location;
    if ~isempty(neighbor_neighbors)
        for ix = 1:numel(neighbor_neighbors)
            neighbor_neighbor_index = neighbor_neighbors(ix);
            [m_mat,n_mat,o_i] = update_neighbor(neighbor_index,neighbor_neighbor_index,0,1,m_mat,n_mat,A,obs_n,o_i);
        end
    end
end