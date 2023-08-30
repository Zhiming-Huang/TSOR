clear;
%% Initialize the network
% load('RWP_input_6_nodes.mat');
communication_Radius = [5,5,5,5,5,5,5,5,5,5];
N = 10 ;
s_input = struct('V_POSITION_X_INTERVAL',[10 50],...%(m)
                 'V_POSITION_Y_INTERVAL',[10 50],...%(m)
                 'V_SPEED_INTERVAL',[2 2],...%(m/s)
                 'V_PAUSE_INTERVAL',[0 1],...%pause time (s)
                 'V_WALK_INTERVAL',[4.00 6.00],...%walk time (s)
                 'V_DIRECTION_INTERVAL',[-180 180],...%(degrees)
                 'SIMULATION_TIME',4300,...%(s)
                 'NB_NODES',8);
s_mobility = Generate_Mobility(s_input);
timeStep = 1;%(s)
% communication_Radius = 8;
[A,A_p,Obs] = generate_Obs_f(s_mobility,s_input,timeStep,communication_Radius);
% [A,A_p,Obs] = generate_Obs_f(s_mobility,s_input,timeStep,communication_Radius);

% upper bound on the node degree
B = max(sum(A,2));
% P(S|i) in hindsight, once the p is known, P(S|i) can be calculated
PSi = cell(N,1);
PSi_p = cell(N,1);
for i = 1:N
    neighbor_S = find(A(i,:) > 0);
    PS = powerset(neighbor_S);
    PSi(i) = {PS};
    p_i = A_p(i,:);
    p_i_c = 1 - p_i;
    p_array = zeros(1,numel(PS));
    for j = 1:numel(PS)
        Z = cell2mat(PS(j));
        Z_c = setdiff(neighbor_S, Z);
        pz = [p_i(Z), p_i_c(Z_c)];
        p_array(j) = prod(pz);
    end
    PSi_p(i) = {p_array};
end
src = 1;
% destination node
dst = N;
% transmission cost
c_i = ones(N,1);
% R
R = 15;
% r
r_i = zeros(N,1);
r_i(N) = R;
% Calculate Vi_star
Vi_star = zeros(N,1);
Vi_star(N) = -R;
Vi_star = calV(Vi_star, A, PSi, PSi_p, N, c_i, r_i);
% Calculate alpha
% alpha = R * sum(1/c_i);
% Calculate delta
% delta = R;
% for i = 1:N
%     neighbor_S = find(A(i,:) > 0);
%     for j = 1:numel(neighbor_S)
%         delta_tmp = abs(Vi_star(i)-Vi_star(neighbor_S(j)));
%         if delta_tmp < delta
%             delta = delta_tmp;
%         end
%     end
% end
% G = N * alpha^2 * B^2 / (2 * delta^2);

%% Other parameter initialized
% number of packets
T = 200;
% generation time of each packet
sigma_t = zeros(T,1);
% arrival rate, according to poisson process
lambda = 0.5;
% EXP_TIMES
EXP_TIMES = 5;

%% Generate packets by poisson process with lambda arrival rate
% T_num = 0;
% n = 0;
% while 1
%     n = n + 1;
%     packet_num = poissrnd(lambda);
%     sigma_t(T_num + 1:T_num+packet_num) = n;
%     T_num = T_num + packet_num;
%     if T_num > T
%         break;
%     end
% end
% sigma_t = sigma_t(1:T);
% n_max = sigma_t(T) + 10;

% while 1
%     n = n + 1;
%     if binornd(1,lambda)
%         sigma_t(T_num) = n;
%         T_num = T_num + 1;
%         if T_num > T
%             break;
%         end
%     end
% end

%% Exp start
Reward_summary = zeros(T,4);
Regret_summary = zeros(T,3);
V2_summary_tsde = zeros(T,N);
V2_summary_dorl = zeros(T,N);
V2_summary_dorl2 = zeros(T,N);
for ext = 1: EXP_TIMES
    disp(ext);
%     Obs = zeros(N,N,n_max);
    % Generate the observations
%     [A,A_p,Obs] = generate_Obs_f(s_mobility,s_input,timeStep,communication_Radius);
%     PSi = cell(N,1);
%     PSi_p = cell(N,1);
%     for i = 1:N
%         neighbor_S = find(A(i,:) > 0);
%         PS = powerset(neighbor_S);
%         PSi(i) = {PS};
%         p_i = A_p(i,:);
%         p_i_c = 1 - p_i;
%         p_array = zeros(1,numel(PS));
%         for j = 1:numel(PS)
%             Z = cell2mat(PS(j));
%             Z_c = setdiff(neighbor_S, Z);
%             pz = [p_i(Z), p_i_c(Z_c)];
%             p_array(j) = prod(pz);
%         end
%         PSi_p(i) = {p_array};
%     end
%     
%     % Calculate Vi_star
%     Vi_star = zeros(N,1);
%     Vi_star(N) = -R;
%     Vi_star = calV(Vi_star, A, PSi, PSi_p, N, c_i, r_i);
    T_num = 0;
    n = 0;
    while 1
        n = n + 1;
        packet_num = poissrnd(lambda);
        sigma_t(T_num + 1:T_num+packet_num) = n;
        T_num = T_num + packet_num;
        if T_num > T
            break;
        end
    end
    sigma_t = sigma_t(1:T);
    n_max = sigma_t(T) + 10;
    
    %% OPT algorithm
    reward_opt = OPT(Vi_star,T,R,N,Obs,sigma_t,n_max,c_i);
    Reward_summary(:,1) = Reward_summary(:,1) + reward_opt;
    
    %% TS_DE algorithm
    [reward_tsde,Vi_hat,PSi_p_hat,v2_tsde] = TSOR_DE(T,R,N,PSi,Obs,sigma_t,n_max,A,A_p,c_i,r_i);
    Reward_summary(:,2) = Reward_summary(:,2) + reward_tsde;
    V2_summary_tsde = V2_summary_tsde + v2_tsde;
    Regret_summary(:,1) = (reward_opt - reward_tsde) + Regret_summary(:,1);
    
    %% DORL with G = 10 algorithm
    [reward_dorl,Vi_hat_dorl,PSi_p_hat_dorl,v2_dorl] = DORL(T,R,N,PSi,Obs,sigma_t,n_max,A,c_i,r_i,10);
    Reward_summary(:,3) = Reward_summary(:,3) + reward_dorl;
    V2_summary_dorl = V2_summary_dorl + v2_dorl;
    Regret_summary(:,2) = (reward_opt - reward_dorl) + Regret_summary(:,2);
    
    %% DORL with G = 50 algorithm
    [reward_dorl2,Vi_hat_dorl2,PSi_p_hat_dorl2,v2_dorl2] = DORL(T,R,N,PSi,Obs,sigma_t,n_max,A,c_i,r_i,50);
    Reward_summary(:,4) = Reward_summary(:,4) + reward_dorl2;
    V2_summary_dorl2 = V2_summary_dorl2 + v2_dorl2;
    Regret_summary(:,3) = (reward_opt - reward_dorl2) + Regret_summary(:,3);
    
    %% TS algorithm
%     [reward_ts,Vi_hat,PSi_p_hat] = TSOR(T,R,N,PSi,Obs,sigma_t,n_max,A,c_i,r_i);
%     Reward_summary_TS = Reward_summary_TS + reward_ts;
%     Regret_summary_TS = (reward_opt - reward_ts) + Regret_summary_TS;
end

% Regret_summary_TS = Regret_summary_TS / EXP_TIMES;
% 
% Packet_average_regret = zeros(T,1);
% for t=1:T
%     Packet_average_regret(t) = sum(Regret_summary_TS(1:t)) / t;
% end
% plot(Packet_average_regret);

Regret_summary = Regret_summary / EXP_TIMES;
Reward_summary = Reward_summary / EXP_TIMES;
V2_summary_tsde = V2_summary_tsde / EXP_TIMES;
V2_summary_dorl = V2_summary_dorl / EXP_TIMES;
V2_summary_dorl2 = V2_summary_dorl2 / EXP_TIMES;
% Packet_average_regret_TSDE = zeros(T,1);
packet_average_regret = zeros(T,3);
packet_average_reward = zeros(T,4);

for t=1:T
    packet_average_regret(t,:) = sum(Regret_summary(1:t,:),1) / t;
    packet_average_reward(t,:) = sum(Reward_summary(1:t,:),1) / t;
end
