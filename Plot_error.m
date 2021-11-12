figure
er=error(:,11:end)-error2(:,11:end);
plot(total_Time(1:T_sim),er)


figure
xxx=xh-xhlin;
plot(total_Time(1:T_sim),xxx(:,1:end-1))

xxx1=[zeros(MatSize,1) xxx];
xxx2=xxx1(1:7,1:end-1)-xxx(1:7,:);

figure
plot(xxx2')