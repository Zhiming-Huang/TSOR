function generate_Obs(s_mobility,s_input,time_step,communication_Radius)
    
    v_t = 0:time_step:s_input.SIMULATION_TIME;
    
    for nodeIndex = 1:s_mobility.NB_NODES
        %Simple interpolation (linear) to get the position, anytime.
        %Remember that "interp1" is the matlab function to use in order to
        %get nodes' position at any continuous time.
        vs_node(nodeIndex).v_x = interp1(s_mobility.VS_NODE(nodeIndex).V_TIME,s_mobility.VS_NODE(nodeIndex).V_POSITION_X,v_t);
        vs_node(nodeIndex).v_y = interp1(s_mobility.VS_NODE(nodeIndex).V_TIME,s_mobility.VS_NODE(nodeIndex).V_POSITION_Y,v_t);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Obs = zeros(s_mobility.NB_NODES,s_mobility.NB_NODES,length(v_t));
    
%     figure;
%     hold on;
%     for nodeIndex = 1:s_mobility.NB_NODES
%         vh_node_pos(nodeIndex) = plot(vs_node(nodeIndex).v_x(1),vs_node(nodeIndex).v_y(1),'*','color',[0.3 0.3 1]);
%     end
%     title(cat(2,'Simulation time (sec): ',num2str(s_mobility.SIMULATION_TIME)));
%     xlabel('X (meters)');
%     ylabel('Y (meters)');
%     title('Radom Waypoint mobility');
%     ht = text(min(vs_node(1).v_x),max(vs_node(1).v_y),cat(2,'Time (sec) = 0'));
%     axis([min(vs_node(1).v_x) max(vs_node(1).v_x) min(vs_node(1).v_y) max(vs_node(1).v_y)]);
%     hold off;
    
    N = s_mobility.NB_NODES;
    contact_sum = zeros(N,N);
    for timeIndex = 1:length(v_t)
%         t = v_t(timeIndex);
%         set(ht,'String',cat(2,'Time (sec) = ',num2str(t,4)));
        obs = zeros(N,N);
        for nodeIndex = 1:s_mobility.NB_NODES
%             set(vh_node_pos(nodeIndex),'XData',vs_node(nodeIndex).v_x(timeIndex),'YData',vs_node(nodeIndex).v_y(timeIndex));
            for neighborIndex = 1:s_mobility.NB_NODES
                if neighborIndex == nodeIndex
                    continue;
                end
                if norm([vs_node(nodeIndex).v_x(timeIndex),vs_node(nodeIndex).v_y(timeIndex)] - [vs_node(neighborIndex).v_x(timeIndex),vs_node(neighborIndex).v_y(timeIndex)]) <= communication_Radius
                    obs(nodeIndex, neighborIndex) = 1;
                    contact_sum(nodeIndex, neighborIndex) = contact_sum(nodeIndex, neighborIndex) + 1;
%                     hold on;
%                     plot([vs_node(nodeIndex).v_x(timeIndex),vs_node(nodeIndex).v_y(timeIndex)],[vs_node(neighborIndex).v_x(timeIndex),vs_node(neighborIndex).v_y(timeIndex)]);
                end
            end
        end
        Obs(:,:,timeIndex) = obs;
%         drawnow;
    end
    
    A_p = contact_sum / length(v_t);
    A = A_p;
    A(A_p > 0) = 1;
end