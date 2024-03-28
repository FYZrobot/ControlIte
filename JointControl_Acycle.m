%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%     控制规划     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%    利用主动关节的正弦振动输入将被动关节调控到期望位置    %%%%%%%%%%%%%%%%%%%%%%%%
close all;clear all; clc;
tic
load CSof7DF;
t0=0.05; %时间步长
n=400; %被动关节规划步数，总步数
wn=1;
sigma=0.9;
kp=wn^2;
kd=2*sigma*wn;

qint=[40,150,60,-60,15,-85,40]'*pi/180;%初始关节角
q=zeros(7,n+1);
q(:,1)=qint;
vq=zeros(7,n+1);
aq=zeros(7,n);
taom=zeros(7,n);
failnum=4;  %故障关节
activenum=3; %主动关节
qd_fail=-110*pi/180;  %被动关节期望角度

controbility=zeros(n,2); %存储运动过程的Mpa的值，观察是否发生奇异
%% 主函数， 利用主动关节的正弦振动输入将被动关节调控到期望位置
omiga=2*pi;
A=zeros(n+1,1);lim_A=0.01;A(1,1)=0.005;
for i=1:n
    q(activenum,i)=q(activenum,1)+A(i,1)*cos(omiga*i*t0);
    vq(activenum,i)=-A(i,1)*omiga*sin(omiga*i*t0);
    D=H(q(:,i));
    C=FCK(q(:,i),vq(:,i));
    controbility(i,1)=D(failnum,activenum)/D(failnum,failnum); %被动关节加速度与主动关节加速度间的关系
    controbility(i,2)=D(failnum,activenum)/(D(failnum,activenum)*D(activenum,failnum)-D(activenum,activenum)*D(failnum,failnum)); %被动关节加速度与主动关节力矩间的关系

    AA=(D(failnum,failnum)*(kp*(qd_fail-q(failnum,i))+kd*(0-vq(failnum,i)))+C(failnum,1))/D(activenum,failnum)/omiga^2*cos(omiga*i*t0); %主动关节的振幅
    %%限定下A的取值范围，防止因在奇异构型附近而计算A值过大（D(activenum,failnum)接近于0值），主动关节正弦输入超出意外之值
%     if abs(AA)<=lim_A
%         A(i+1,1)=AA;
%     else
%         A(i+1,1)=lim_A*sign(AA);
%     end
    A(i+1,1)=AA;
    aq(activenum,i)=-A(i+1,1)*omiga^2*cos(omiga*i*t0); % 主动关节加速度
    
    aq(failnum,i)=(-aq(activenum,i)*D(failnum,activenum)-C(failnum,1))/D(failnum,failnum); % 被动关节加速度
    vq(failnum,i+1)=vq(failnum,i)+aq(failnum,i)*t0;
    q(failnum,i+1)=q(failnum,i)+vq(failnum,i)*t0;
    taom(:,i)=H(q(:,i))*aq(:,i)+FCK(q(:,i),vq(:,i));
end

%% 关节曲线(简) %%%%%%
tt=1:n;
%%%期望角度
desire_qp=zeros(n,1);
for i=1:n
    desire_qp(i,1)=qd_fail;
end
figure(1) %角度
plot(tt*t0,q(activenum,tt)*180/pi,'k-',tt*t0,q(failnum,tt)*180/pi,'k--',tt*t0,desire_qp(tt,1)*180/pi,'k:','linewidth',2);
xlabel('Time (s)');
ylabel('Angle (°)');
legend('The Active Joint','The Passive Joint','The Desire Angle');
title('角度');
hold on;
grid on;
% figure(1) %角度
% plot(tt*t0,q(failnum,tt)*180/pi,'k-',tt*t0,desire_qp(tt,1)*180/pi,'k--','linewidth',2);
% xlabel('Time (s)');
% ylabel('Angle (°)');
% legend('The free-swinging joint','The desire angle');
% hold on;
% grid on;
figure(2) %角速度
plot(tt*t0,vq(activenum,tt)*180/pi,'k-',tt*t0,vq(failnum,tt)*180/pi,'k--','linewidth',2);
xlabel('Time (s)');
ylabel('Vlocity (°/s)');
legend('The Active Joint','The Passive Joint');
title('角速度');
hold on;
grid on;
figure(3) %角加速度
plot(tt*t0,aq(activenum,tt)*180/pi,'k-.',tt*t0,aq(failnum,tt)*180/pi,'k-','linewidth',2);
xlabel('Time (s)');
ylabel('Angular acceleration (°/s^2)');
legend('The actuated joint','The free-swinging joint');
% title('角加速度');
hold on;
grid on;
figure(4) %关节力矩
plot(tt*t0,taom(activenum,tt),'k-',tt*t0,taom(failnum,tt),'k--','linewidth',2);
xlabel('Time (s)');
ylabel('Torque (NM)');
legend('The Active Joint','The Passive Joint');
title('力矩');
hold on;
grid on;
figure(5) %耦合程度
plot(tt*t0,controbility(tt,1),'k-',tt*t0,controbility(tt,2),'k--','linewidth',2);
xlabel('Time (s)');
legend('加速度关系','加速度与力矩间关系');
title('耦合程度');
hold on;
grid on;
figure(6) %A
plot(tt*t0,A(tt,1),'k-','linewidth',2);
xlabel('Time (s)');
ylabel('Amplitude, A ');
% title('A');
hold on;
grid on;

toc