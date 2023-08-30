marker_size = 10;
font_size = 20;
fig_num = 0;
sce_no = '1';

fig_num = fig_num + 1;
fig = figure(fig_num);
plot_T_num = T;
plot_T = 1:plot_T_num;
plot(packet_average_regret(:,1),'Marker','h', 'MarkerSize',marker_size,'MarkerIndices',1:floor(plot_T_num/10):plot_T_num);
hold on;
plot(packet_average_regret(:,2),'Marker','o', 'MarkerSize',marker_size,'MarkerIndices',1:floor(plot_T_num/10):plot_T_num);
hold on;
plot(packet_average_regret(:,3),'Marker','*', 'MarkerSize',marker_size,'MarkerIndices',1:floor(plot_T_num/10):plot_T_num);
hold on;
plot(zeros(T,1),'Marker','>', 'MarkerSize',marker_size,'MarkerIndices',1:floor(plot_T_num/10):plot_T_num,'color',[1.00,0.07,0.65]);
hold on;
legend('TSOR','DORL, G=10','DORL, G=50', 'OPT')
grid on;
xlabel('# of Packets','FontSize',font_size);
ylabel('Packet-averaged regret','FontSize',font_size);
set(gca, 'FontSize', font_size);

set(fig, 'PaperUnits', 'inches');
set(fig, 'PaperSize', [8.5 6.25]);
set(fig, 'PaperPositionMode', 'manual');
set(fig, 'PaperPosition', [0 0 8.5 6.25]);
set(fig, 'renderer', 'painters');
print(fig, '-depsc2', ['result_sce_',sce_no,'_regret.eps']);

fig_num = fig_num + 1;
fig = figure(fig_num);
plot_T_num = T;
plot_T = 1:plot_T_num;
plot(packet_average_reward(:,2),'Marker','h', 'MarkerSize',marker_size,'MarkerIndices',1:floor(plot_T_num/10):plot_T_num);
hold on;
plot(packet_average_reward(:,3),'Marker','o', 'MarkerSize',marker_size,'MarkerIndices',1:floor(plot_T_num/10):plot_T_num);
hold on;
plot(packet_average_reward(:,4),'Marker','*', 'MarkerSize',marker_size,'MarkerIndices',1:floor(plot_T_num/10):plot_T_num);
hold on;
plot(packet_average_reward(:,1),'Marker','>', 'MarkerSize',marker_size,'MarkerIndices',1:floor(plot_T_num/10):plot_T_num,'color',[1.00,0.07,0.65]);
hold on;
legend({'TSOR','DORL, G=10','DORL, G=50','OPT'},'Location', 'northeast');
grid on;
xlabel('# of Packets','FontSize',font_size);
ylabel('Packet-averaged reward','FontSize',font_size);
set(gca, 'FontSize', font_size);

set(fig, 'PaperUnits', 'inches');
set(fig, 'PaperSize', [8.5 6.25]);
set(fig, 'PaperPositionMode', 'manual');
set(fig, 'PaperPosition', [0 0 8.5 6.25]);
set(fig, 'renderer', 'painters');
print(fig, '-depsc2', ['result_sce_',sce_no,'_reward.eps']);

fig_num = fig_num + 1;
fig = figure(fig_num);
plot_T_num = T;
plot_T = 1:plot_T_num;
plot(R-packet_average_reward(:,2),'Marker','h', 'MarkerSize',marker_size,'MarkerIndices',1:floor(plot_T_num/10):plot_T_num);
hold on;
plot(R-packet_average_reward(:,3),'Marker','o', 'MarkerSize',marker_size,'MarkerIndices',1:floor(plot_T_num/10):plot_T_num);
hold on;
plot(R-packet_average_reward(:,4),'Marker','*', 'MarkerSize',marker_size,'MarkerIndices',1:floor(plot_T_num/10):plot_T_num);
hold on;
plot(R-packet_average_reward(:,1),'Marker','>', 'MarkerSize',marker_size,'MarkerIndices',1:floor(plot_T_num/10):plot_T_num,'color',[1.00,0.07,0.65]);
hold on;
legend({'TSOR','DORL, G=10','DORL, G=50','OPT'},'Location', 'southeast');
grid on;
xlabel('# of Packets','FontSize',font_size);
ylabel('Packet-averaged cost','FontSize',font_size);
set(gca, 'FontSize', font_size);

set(fig, 'PaperUnits', 'inches');
set(fig, 'PaperSize', [8.5 6.25]);
set(fig, 'PaperPositionMode', 'manual');
set(fig, 'PaperPosition', [0 0 8.5 6.25]);
set(fig, 'renderer', 'painters');
print(fig, '-depsc2', ['result_sce_',sce_no,'_cost.eps']);

fig_num = fig_num + 1;
fig = figure(fig_num);
plot_T_num = T;
plot_T = 1:plot_T_num;
plot(V2_summary_tsde(:,1),'Marker','h', 'MarkerSize',marker_size,'MarkerIndices',1:floor(plot_T_num/10):plot_T_num);
hold on;
plot(V2_summary_dorl(:,1),'Marker','o', 'MarkerSize',marker_size,'MarkerIndices',1:floor(plot_T_num/10):plot_T_num);
hold on;
plot(V2_summary_dorl2(:,1),'Marker','*', 'MarkerSize',marker_size,'MarkerIndices',1:floor(plot_T_num/10):plot_T_num);
hold on;
plot(ones(plot_T_num,1)*Vi_star(1),'Marker','>', 'MarkerSize',marker_size,'MarkerIndices',1:floor(plot_T_num/10):plot_T_num,'color',[1.00,0.07,0.65]);
hold on;
legend({'TSOR','DORL, G=10','DORL, G=50','OPT'},'Location', 'best');
grid on;
xlabel('# of Packets','FontSize',font_size);
ylabel('Estimated value of source node','FontSize',font_size);
set(gca, 'FontSize', font_size);

set(fig, 'PaperUnits', 'inches');
set(fig, 'PaperSize', [8.5 6.25]);
set(fig, 'PaperPositionMode', 'manual');
set(fig, 'PaperPosition', [0 0 8.5 6.25]);
set(fig, 'renderer', 'painters');
print(fig, '-depsc2', ['result_sce_',sce_no,'_v2_summary_n1.eps']);

fig_num = fig_num + 1;
fig = figure(fig_num);
plot_T_num = T;
plot_T = 1:plot_T_num;
plot(V2_summary_tsde(:,2),'Marker','h', 'MarkerSize',marker_size,'MarkerIndices',1:floor(plot_T_num/10):plot_T_num);
hold on;
plot(V2_summary_dorl(:,2),'Marker','o', 'MarkerSize',marker_size,'MarkerIndices',1:floor(plot_T_num/10):plot_T_num);
hold on;
plot(V2_summary_dorl2(:,2),'Marker','*', 'MarkerSize',marker_size,'MarkerIndices',1:floor(plot_T_num/10):plot_T_num);
hold on;
plot(ones(plot_T_num,1)*Vi_star(2),'Marker','>', 'MarkerSize',marker_size,'MarkerIndices',1:floor(plot_T_num/10):plot_T_num,'color',[1.00,0.07,0.65]);
hold on;
legend({'TSOR','DORL, G=10','DORL, G=50','OPT'},'Location', 'best');
grid on;
xlabel('# of Packets','FontSize',font_size);
ylabel('Estimated value of node 2','FontSize',font_size);
set(gca, 'FontSize', font_size);

set(fig, 'PaperUnits', 'inches');
set(fig, 'PaperSize', [8.5 6.25]);
set(fig, 'PaperPositionMode', 'manual');
set(fig, 'PaperPosition', [0 0 8.5 6.25]);
set(fig, 'renderer', 'painters');
print(fig, '-depsc2', ['result_sce_',sce_no,'_v2_summary_n2.eps']);

% % 
% % % fig_num = fig_num + 1;
% % % fig = figure(fig_num);
% % % plot_T_num = T;
% % % plot_T = 1:plot_T_num;
% % % for i = 1:N-1
% % %     plot(V2_summary_tsde(:,i),'Marker','h', 'linestyle', 'None', 'MarkerSize',marker_size,'MarkerIndices',1:floor(plot_T_num/10):plot_T_num);
% % %     hold on;
% % %     plot(ones(plot_T_num,1)*Vi_star(i));
% % %     hold on;
% % % end
% % % legend({'TSOR','DORL, G=10','DORL, G=50','Target'},'Location', 'northeast');
% % % grid on;
% % % xlabel('# of Packets','FontSize',font_size);
% % % ylabel('Estimated value','FontSize',font_size);
% % % set(gca, 'FontSize', font_size);
