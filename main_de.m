clear;
%% Initialize the network
% number of nodes
% N = 4;
% s = [1,1,2,2,3];
% t = [2,3,3,4,4];
% p = [0.6,0.4,0.5,0.3,0.6];
N = 6;
s = [1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,5,5,5,5];
t = [2,3,4,5,1,3,4,5,6,1,2,4,5,6,1,2,3,5,1,2,3,4];
%%% link probability
load('/Users/imote/Dropbox/Mobihoc''20/MATLAB code/p_6.mat');
% p = ones(1,numel(s)) * 0.1 .* randi(9,1,numel(s)); source node
src = 1;
% destination node
dst = N;
% transmission cost
c_i = ones(N,1);
% R
R = 10;
% r
r_i = zeros(N,1);
r_i(N) = R;
% s = [1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5];
% t = [1,2,3,4,5,1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,1,2,3,4,5];
Gr = digraph(s,t,p);
% adjacency matrix
A = full(adjacency(Gr));
% adjacency matrix with weight
A_p = full(adjacency(Gr,'weighted'));
% plot(Gr,'Layout','force','EdgeLabel',Gr.Edges.Weight);
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
% Calculate Vi_star
Vi_star = zeros(N,1);
Vi_star(N) = -R;
Vi_star = calV(Vi_star, A, PSi, PSi_p, N, c_i, r_i);
% Calculate alpha
alpha = R * sum(1/c_i);
% Calculate delta
delta = R;
for i = 1:N
    neighbor_S = find(A(i,:) > 0);
    for j = 1:numel(neighbor_S)
        delta_tmp = abs(Vi_star(i)-Vi_star(neighbor_S(j)));
        if delta_tmp < delta
            delta = delta_tmp;
        end
    end
end
% G = N * alpha^2 * B^2 / (2 * delta^2);

%% Other parameter initialized
% number of packets
T = 2000;
% generation time of each packet
sigma_t = zeros(T,1);
% arrival rate, according to poisson process
lambda = 0.5;
% EXP_TIMES
EXP_TIMES = 1;

%% Generate packets by poisson process with lambda arrival rate
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
n_max = sigma_t(T) + 10;
Reward_summary_TS = zeros(T,1);
Reward_summary_OPT = zeros(T,1);
Regret_summary_TS = zeros(T,1);
for ext = 1: EXP_TIMES
    disp(ext);
    Obs = zeros(N,N,n_max);
    % Generate the observations
    for i = 1:n_max
        Obs(:,:,i) = binornd(1,A_p);
    end
    %% OPT algorithm
    reward_opt = OPT(Vi_star,T,R,N,Obs,sigma_t,n_max,c_i);
    Reward_summary_OPT = Reward_summary_OPT + reward_opt;
    
    %% TS algorithm
    [reward_ts,Vi_hat] = TSOR(T,R,N,PSi,Obs,sigma_t,n_max,A,c_i,r_i);
    Reward_summary_TS = Reward_summary_TS + reward_ts;
    Regret_summary_TS = (reward_opt - reward_ts) + Regret_summary_TS;
end

Regret_summary_TS = Regret_summary_TS / EXP_TIMES;

% Packet_average_cost = zeros(T,1);
% for t=1:T
%     Packet_average_cost(t) = sum(Cost_summary_TS(1:t)) / t;
% end
% plot(Packet_average_cost);
Packet_average_regret = zeros(T,1);
for t=1:T
    Packet_average_regret(t) = sum(Regret_summary_TS(1:t)) / t;
end
plot(Packet_average_regret);
