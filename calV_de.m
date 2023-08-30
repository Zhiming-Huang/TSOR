function Vi = calV_de(Vi, A, PSi, PSi_p, c_i, r_i,Set_A_flag,Set_X_flag)

% Set_A_flag = zeros(N,1);
% Set_A_flag(N) = 1;
Set_A = find(Set_A_flag == 1);
Set_X = find(Set_X_flag == 1);
while numel(Set_X) > 0
    m = numel(Set_X);
    for i = 1:m
%         find(A(Set_X(i),:) > 0)
%         if isempty(intersect(find(A(Set_X(i),:)>0),Set_A))
%             continue;
%         end
        if all(ismember(find(A(Set_X(i),:)>0),Set_A) == 0)
            continue;
        end
        PSi_i = PSi(Set_X(i));
        PSi_i = PSi_i{1,1};
        PSi_i_p = PSi_p(Set_X(i));
        PSi_i_p = PSi_i_p{1,1};
        temp_sum_up = 0;
        temp_sum_down = 0;
        for j = 2:numel(PSi_i)
            psi = cell2mat(PSi_i(j));
            if all(ismember(psi,Set_A) == 0)
                continue;
            end
            min_V = 0;
            for k = 1:numel(psi)
                if Vi(psi(k)) < min_V
                    min_V = Vi(psi(k));
                end
            end
            temp_sum_up = PSi_i_p(j) * min_V + temp_sum_up;
            temp_sum_down = PSi_i_p(j) + temp_sum_down;
        end
        min_right = (c_i(Set_X(i)) + temp_sum_up) / temp_sum_down;
        Vi(Set_X(i)) = min([-r_i(Set_X(i)),min_right]);
    end
    [min_value,min_index] = min(Vi(Set_X));
    Set_A_flag(Set_X(min_index)) = 1;
    Set_A = find(Set_A_flag == 1);
    Set_X_flag(Set_X(min_index)) = 0;
    Set_X = find(Set_X_flag == 1);
%     Set_X = setdiff(1:N, Set_A);
end