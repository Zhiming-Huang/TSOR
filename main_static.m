clear;
%% Initialize the network
% number of nodes
% N = 4;
% s = [1,1,2,2,3,2,3,3,4,4];
% t = [2,3,1,3,2,4,1,4,2,3];
% p = [0.6,0.4,0.6,0.5,0.5,0.3,0.4,0.6,0.5,0.5];
N = 6;
s = [1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,5,5,5,5,6,6];
t = [2,3,4,5,1,3,4,5,6,1,2,4,5,6,1,2,3,5,1,2,3,4,2,3];
% %%% link probability
% load('p_6d.mat');
% p = ones(1,numel(s)) * 0.1 .* randi(9,1,numel(s));
% N = 16;
% s = [1,1,2,2,2,3,3,3,4,4,5,5,5,6,6,6,6,7,7,7,7,8,8,8,9,9,9,10,10,10,10,11,11,11,11,12,12,12,13,13,14,14,14,15,15,15,16,16];
% t = [2,5,1,3,6,2,4,7,3,8,1,6,9,2,5,7,10,3,6,8,11,4,7,12,5,10,13,6,9,11,14,7,10,12,15,8,11,16,9,14,10,13,15,11,14,16,12,15];
p = ones(1,numel(s)) * 0.1 .* randi(9,1,numel(s));
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
% plotgraph = plot(Gr,'Layout','force','EdgeLabel',Gr.Edges.Weight);
Gr_undirected = graph(A);
% plotgraph = plot(Gr_undirected);
% plotgraph.XData(src) = 0;
% plotgraph.YData(src) = 0.5;
% plotgraph.XData(dst) = 3;
% plotgraph.YData(dst) = 0.5;
% plotgraph.XData(2:3) = 2;
% plotgraph.YData(2:3) = linspace(0,0.4,2);
% plotgraph.XData(4:5) = 1;
% plotgraph.YData(4:5) = linspace(0.6,1,2);
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
T = 3000;
% generation time of each packet
sigma_t = zeros(T,1);
% arrival rate, according to poisson process
lambda = 0.5;
% EXP_TIMES
EXP_TIMES = 100;

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
Reward_summary = zeros(T,4);
Regret_summary = zeros(T,3);
V2_summary_tsde = zeros(T,N);
V2_summary_dorl = zeros(T,N);
V2_summary_dorl2 = zeros(T,N);
for ext = 1: EXP_TIMES
    disp(ext);
    Obs = zeros(N,N,n_max);
    % Generate the observations
    for i = 1:n_max
        Obs(:,:,i) = binornd(1,A_p);
    end
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

% plot(packet_average_regret);

