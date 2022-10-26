% This script returns figures plotting the input, reference, and obtained 
% output trajectories for the given DASH specs.

function DASH_plot_output(i, x, v, g_b, g_a, maxtrajsize, p)

ctext = [];
switch p.io
    case 1
        ctext = [ctext 'Simple_Err_'];
    case 2
        ctext = [ctext 'Complex_Err_'];
end
switch p.c
    case 0
        ctext = [ctext 'Mount_Swing_Lift'];
    case 1
        ctext = [ctext 'Mount_Lift'];
    case 2
        ctext = [ctext 'Mount_Swing'];
end

circ1 = (1/255)*[215, 25, 28];
circ3 = (1/255)*[44, 123, 182];
% circ2 = (1/255)*[96, 0, 220];
% circ4 = (1/255)*[35,139,69];
out = DASH_animate(x', v, p);
h_i_t = out.h_i_t;
h_o_t = out.h_o_t;
t = out.t;
i_mean_x = mean(h_i_t(1,4,:)); i_mean_y = mean(h_i_t(2,4,:)); i_mean_z = mean(h_i_t(3,4,:));
o_mean_x = mean(h_o_t(1,4,:)); o_mean_y = mean(h_o_t(2,4,:)); o_mean_z = mean(h_o_t(3,4,:));
xi = reshape(h_i_t(1,4,:),1,length(t))-i_mean_x;
yi = reshape(h_i_t(2,4,:),1,length(t))-i_mean_y;
zi = reshape(h_i_t(3,4,:),1,length(t))-i_mean_z;
xo = reshape(h_o_t(1,4,:),1,length(t))-o_mean_x;
yo = reshape(h_o_t(2,4,:),1,length(t))-o_mean_y;
zo = reshape(h_o_t(3,4,:),1,length(t))-o_mean_z;
r_i = 0.25*(max(yi) - min(yi) + max(zi) - min(zi)); 
xo_r = zeros(size(t));
yo_r = g_b*r_i*sin(2*pi*t - pi/2);
zo_r = g_a*r_i*sin(2*pi*t - pi);
figure('units','pixels','position',[0 0 1920 1080],'Color','w')
subplot(1,3,1)
plot3(xi(1),zi(1),yi(1),'o','LineWidth',1.2,'Color',circ1);
hold on; grid on; axis equal;
zlim([-maxtrajsize maxtrajsize]); ylim([-maxtrajsize maxtrajsize]);
plot3(xi,zi,yi,'Color',circ1,'LineWidth',1.2);
view(-90,0);
plot3(xo(1),zo(1),yo(1),'o','LineWidth',1.2,'Color',circ3);
plot3(xo,zo,yo,'Color',circ3,'LineWidth',1.2);
plot3(xo_r(1),zo_r(1),yo_r(1),'o--','LineWidth',1.0,'Color',circ3);
plot3(xo_r,zo_r,yo_r,'--','Color',circ3,'LineWidth',1.0);
xlabel('x');
ylabel('y');
zlabel('z');
subplot(1,3,2)
plot3(xi(1),zi(1),yi(1),'o','LineWidth',1.2,'Color',circ1);
hold on; grid on; axis equal;
zlim([-maxtrajsize maxtrajsize]); ylim([-maxtrajsize maxtrajsize]);
view([0 0]);
plot3(xi,zi,yi,'Color',circ1,'LineWidth',1.2);
plot3(xo(1),zo(1),yo(1),'o','LineWidth',1.2,'Color',circ3);
plot3(xo,zo,yo,'Color',circ3,'LineWidth',1.2);
plot3(xo_r(1),zo_r(1),yo_r(1),'o--','LineWidth',1.0,'Color',circ3);
plot3(xo_r,zo_r,yo_r,'--','Color',circ3,'LineWidth',1.0);
xlabel('x');
ylabel('y');
zlabel('z');
subplot(1,3,3)
plot3(xi(1),zi(1),yi(1),'o','LineWidth',1.2,'Color',circ1);
hold on; grid on; axis equal;
zlim([-maxtrajsize maxtrajsize]); ylim([-maxtrajsize maxtrajsize]);
view([-20 15]);
plot3(xi,zi,yi,'Color',circ1,'LineWidth',1.2);
plot3(xo(1),zo(1),yo(1),'o','LineWidth',1.2,'Color',circ3);
plot3(xo,zo,yo,'Color',circ3,'LineWidth',1.2);
plot3(xo_r(1),zo_r(1),yo_r(1),'o--','LineWidth',1.0,'Color',circ3);
plot3(xo_r,zo_r,yo_r,'--','Color',circ3,'LineWidth',1.0);
xlabel('x');
ylabel('y');
zlabel('z');
sgtitle(['case ' num2str(i) ctext]);
% isoF = 1;
% out_lift = DASH_animate(x', v, isoF);
% h_i_t = out_lift.h_i_t;
% h_o_t = out_lift.h_o_t;
% t = out_lift.t;
% i_mean_x = mean(h_i_t(1,4,:)); i_mean_y = mean(h_i_t(2,4,:)); i_mean_z = mean(h_i_t(3,4,:));
% o_mean_x = mean(h_o_t(1,4,:)); o_mean_y = mean(h_o_t(2,4,:)); o_mean_z = mean(h_o_t(3,4,:));
% xi = reshape(h_i_t(1,4,:),1,length(t))-i_mean_x;
% yi = reshape(h_i_t(2,4,:),1,length(t))-i_mean_y;
% zi = reshape(h_i_t(3,4,:),1,length(t))-i_mean_z;
% xo = reshape(h_o_t(1,4,:),1,length(t))-o_mean_x;
% yo = reshape(h_o_t(2,4,:),1,length(t))-o_mean_y;
% zo = reshape(h_o_t(3,4,:),1,length(t))-o_mean_z;
% r_i = 0.25*(max(yi) - min(yi) + max(zi) - min(zi)); 
% xo_r = zeros(size(t));
% yo_r = g_b*r_i*sin(2*pi*t - pi/2);
% zo_r = g_a*r_i*sin(2*pi*t - pi);
% figure('units','pixels','position',[0 0 1920 1080],'Color','w')
% subplot(2,3,1)
% plot3(xi(1),zi(1),yi(1),'o','LineWidth',1.2,'Color',circ1);
% hold on; grid on; axis equal;
% % zlim([-maxtrajsize maxtrajsize]); ylim([-maxtrajsize maxtrajsize]);
% plot3(xi,zi,yi,'Color',circ1,'LineWidth',1.2);
% view(-90,0);
% plot3(xo(1),zo(1),yo(1),'o','LineWidth',1.2,'Color',circ3);
% plot3(xo,zo,yo,'Color',circ3,'LineWidth',1.2);
% plot3(xo_r(1),zo_r(1),yo_r(1),'o--','LineWidth',1.0,'Color',circ3);
% plot3(xo_r,zo_r,yo_r,'--','Color',circ3,'LineWidth',1.0);
% xlabel('x');
% ylabel('y');
% zlabel('z');
% subplot(2,3,2)
% plot3(xi(1),zi(1),yi(1),'o','LineWidth',1.2,'Color',circ1);
% hold on; grid on; axis equal;
% % zlim([-maxtrajsize maxtrajsize]); ylim([-maxtrajsize maxtrajsize]);
% view([0 0]);
% plot3(xi,zi,yi,'Color',circ1,'LineWidth',1.2);
% plot3(xo(1),zo(1),yo(1),'o','LineWidth',1.2,'Color',circ3);
% plot3(xo,zo,yo,'Color',circ3,'LineWidth',1.2);
% plot3(xo_r(1),zo_r(1),yo_r(1),'o--','LineWidth',1.0,'Color',circ3);
% plot3(xo_r,zo_r,yo_r,'--','Color',circ3,'LineWidth',1.0);
% xlabel('x');
% ylabel('y');
% zlabel('z');
% subplot(2,3,3)
% plot3(xi(1),zi(1),yi(1),'o','LineWidth',1.2,'Color',circ1);
% hold on; grid on; axis equal;
% % zlim([-maxtrajsize maxtrajsize]); ylim([-maxtrajsize maxtrajsize]);
% view([-20 15]);
% plot3(xi,zi,yi,'Color',circ1,'LineWidth',1.2);
% plot3(xo(1),zo(1),yo(1),'o','LineWidth',1.2,'Color',circ3);
% plot3(xo,zo,yo,'Color',circ3,'LineWidth',1.2);
% plot3(xo_r(1),zo_r(1),yo_r(1),'o--','LineWidth',1.0,'Color',circ3);
% plot3(xo_r,zo_r,yo_r,'--','Color',circ3,'LineWidth',1.0);
% xlabel('x');
% ylabel('y');
% zlabel('z');
% isoF = 2;
% out_swing = DASH_animate(x', v, isoF);
% h_i_t = out_swing.h_i_t;
% h_o_t = out_swing.h_o_t;
% t = out_swing.t;
% i_mean_x = mean(h_i_t(1,4,:)); i_mean_y = mean(h_i_t(2,4,:)); i_mean_z = mean(h_i_t(3,4,:));
% o_mean_x = mean(h_o_t(1,4,:)); o_mean_y = mean(h_o_t(2,4,:)); o_mean_z = mean(h_o_t(3,4,:));
% xi = reshape(h_i_t(1,4,:),1,length(t))-i_mean_x;
% yi = reshape(h_i_t(2,4,:),1,length(t))-i_mean_y;
% zi = reshape(h_i_t(3,4,:),1,length(t))-i_mean_z;
% xo = reshape(h_o_t(1,4,:),1,length(t))-o_mean_x;
% yo = reshape(h_o_t(2,4,:),1,length(t))-o_mean_y;
% zo = reshape(h_o_t(3,4,:),1,length(t))-o_mean_z;
% r_i = 0.25*(max(yi) - min(yi) + max(zi) - min(zi)); 
% xo_r = zeros(size(t));
% yo_r = g_b*r_i*sin(2*pi*t - pi/2);
% zo_r = g_a*r_i*sin(2*pi*t - pi);
% subplot(2,3,4)
% plot3(xi(1),zi(1),yi(1),'o','LineWidth',1.2,'Color',circ1);
% hold on; grid on; axis equal;
% % zlim([-maxtrajsize maxtrajsize]); ylim([-maxtrajsize maxtrajsize]);
% plot3(xi,zi,yi,'Color',circ1,'LineWidth',1.2);
% view(-90,0);
% plot3(xo(1),zo(1),yo(1),'o','LineWidth',1.2,'Color',circ3);
% plot3(xo,zo,yo,'Color',circ3,'LineWidth',1.2);
% plot3(xo_r(1),zo_r(1),yo_r(1),'o--','LineWidth',1.0,'Color',circ3);
% plot3(xo_r,zo_r,yo_r,'--','Color',circ3,'LineWidth',1.0);
% subplot(2,3,5)
% plot3(xi(1),zi(1),yi(1),'o','LineWidth',1.2,'Color',circ1);
% hold on; grid on; axis equal;
% % zlim([-maxtrajsize maxtrajsize]); ylim([-maxtrajsize maxtrajsize]);
% view([0 0]);
% plot3(xi,zi,yi,'Color',circ1,'LineWidth',1.2);
% plot3(xo(1),zo(1),yo(1),'o','LineWidth',1.2,'Color',circ3);
% plot3(xo,zo,yo,'Color',circ3,'LineWidth',1.2);
% plot3(xo_r(1),zo_r(1),yo_r(1),'o--','LineWidth',1.0,'Color',circ3);
% plot3(xo_r,zo_r,yo_r,'--','Color',circ3,'LineWidth',1.0);
% subplot(2,3,6)
% plot3(xi(1),zi(1),yi(1),'o','LineWidth',1.2,'Color',circ1);
% hold on; grid on; axis equal;
% % zlim([-maxtrajsize maxtrajsize]); ylim([-maxtrajsize maxtrajsize]);
% view([-20 15]);
% plot3(xi,zi,yi,'Color',circ1,'LineWidth',1.2);
% plot3(xo(1),zo(1),yo(1),'o','LineWidth',1.2,'Color',circ3);
% plot3(xo,zo,yo,'Color',circ3,'LineWidth',1.2);
% plot3(xo_r(1),zo_r(1),yo_r(1),'o--','LineWidth',1.0,'Color',circ3);
% plot3(xo_r,zo_r,yo_r,'--','Color',circ3,'LineWidth',1.0);
% sgtitle(['case '  num2str(i) ' (row 1: lift alone; row 2: swing alone)']);

end