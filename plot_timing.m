
figure
plot(total_Time(1:T_sim),status)
title(''),
xlabel('Time [s]'),ylabel('Temperature $^\circ C$')
figp

figure
plot(total_Time(1:T_sim),t_quad)
hold on 
plot(total_Time(1:T_sim),t_mpcqp)
title('Computational time per iteration'),
xlabel('Time [s]'),ylabel('Computational time [s]')
figp
legend('Quad','mpcqp')
ylim([0 0.1])
t_total_quad=sum(t_quad)
t_total_mpcqp=sum(t_mpcqp)