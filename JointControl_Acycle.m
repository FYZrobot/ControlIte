%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%     ���ƹ滮     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%    ���������ؽڵ����������뽫�����ؽڵ��ص�����λ��    %%%%%%%%%%%%%%%%%%%%%%%%
close all;clear all; clc;
tic
load CSof7DF;
t0=0.05; %ʱ�䲽��
n=400; %�����ؽڹ滮�������ܲ���
wn=1;
sigma=0.9;
kp=wn^2;
kd=2*sigma*wn;

qint=[40,150,60,-60,15,-85,40]'*pi/180;%��ʼ�ؽڽ�
q=zeros(7,n+1);
q(:,1)=qint;
vq=zeros(7,n+1);
aq=zeros(7,n);
taom=zeros(7,n);
failnum=4;  %���Ϲؽ�
activenum=3; %�����ؽ�
qd_fail=-110*pi/180;  %�����ؽ������Ƕ�

controbility=zeros(n,2); %�洢�˶����̵�Mpa��ֵ���۲��Ƿ�������
%% �������� ���������ؽڵ����������뽫�����ؽڵ��ص�����λ��
omiga=2*pi;
A=zeros(n+1,1);lim_A=0.01;A(1,1)=0.005;
for i=1:n
    q(activenum,i)=q(activenum,1)+A(i,1)*cos(omiga*i*t0);
    vq(activenum,i)=-A(i,1)*omiga*sin(omiga*i*t0);
    D=H(q(:,i));
    C=FCK(q(:,i),vq(:,i));
    controbility(i,1)=D(failnum,activenum)/D(failnum,failnum); %�����ؽڼ��ٶ��������ؽڼ��ٶȼ�Ĺ�ϵ
    controbility(i,2)=D(failnum,activenum)/(D(failnum,activenum)*D(activenum,failnum)-D(activenum,activenum)*D(failnum,failnum)); %�����ؽڼ��ٶ��������ؽ����ؼ�Ĺ�ϵ

    AA=(D(failnum,failnum)*(kp*(qd_fail-q(failnum,i))+kd*(0-vq(failnum,i)))+C(failnum,1))/D(activenum,failnum)/omiga^2*cos(omiga*i*t0); %�����ؽڵ����
    %%�޶���A��ȡֵ��Χ����ֹ�������칹�͸���������Aֵ����D(activenum,failnum)�ӽ���0ֵ���������ؽ��������볬������ֵ֮
%     if abs(AA)<=lim_A
%         A(i+1,1)=AA;
%     else
%         A(i+1,1)=lim_A*sign(AA);
%     end
    A(i+1,1)=AA;
    aq(activenum,i)=-A(i+1,1)*omiga^2*cos(omiga*i*t0); % �����ؽڼ��ٶ�
    
    aq(failnum,i)=(-aq(activenum,i)*D(failnum,activenum)-C(failnum,1))/D(failnum,failnum); % �����ؽڼ��ٶ�
    vq(failnum,i+1)=vq(failnum,i)+aq(failnum,i)*t0;
    q(failnum,i+1)=q(failnum,i)+vq(failnum,i)*t0;
    taom(:,i)=H(q(:,i))*aq(:,i)+FCK(q(:,i),vq(:,i));
end

%% �ؽ�����(��) %%%%%%
tt=1:n;
%%%�����Ƕ�
desire_qp=zeros(n,1);
for i=1:n
    desire_qp(i,1)=qd_fail;
end
figure(1) %�Ƕ�
plot(tt*t0,q(activenum,tt)*180/pi,'k-',tt*t0,q(failnum,tt)*180/pi,'k--',tt*t0,desire_qp(tt,1)*180/pi,'k:','linewidth',2);
xlabel('Time (s)');
ylabel('Angle (��)');
legend('The Active Joint','The Passive Joint','The Desire Angle');
title('�Ƕ�');
hold on;
grid on;
% figure(1) %�Ƕ�
% plot(tt*t0,q(failnum,tt)*180/pi,'k-',tt*t0,desire_qp(tt,1)*180/pi,'k--','linewidth',2);
% xlabel('Time (s)');
% ylabel('Angle (��)');
% legend('The free-swinging joint','The desire angle');
% hold on;
% grid on;
figure(2) %���ٶ�
plot(tt*t0,vq(activenum,tt)*180/pi,'k-',tt*t0,vq(failnum,tt)*180/pi,'k--','linewidth',2);
xlabel('Time (s)');
ylabel('Vlocity (��/s)');
legend('The Active Joint','The Passive Joint');
title('���ٶ�');
hold on;
grid on;
figure(3) %�Ǽ��ٶ�
plot(tt*t0,aq(activenum,tt)*180/pi,'k-.',tt*t0,aq(failnum,tt)*180/pi,'k-','linewidth',2);
xlabel('Time (s)');
ylabel('Angular acceleration (��/s^2)');
legend('The actuated joint','The free-swinging joint');
% title('�Ǽ��ٶ�');
hold on;
grid on;
figure(4) %�ؽ�����
plot(tt*t0,taom(activenum,tt),'k-',tt*t0,taom(failnum,tt),'k--','linewidth',2);
xlabel('Time (s)');
ylabel('Torque (NM)');
legend('The Active Joint','The Passive Joint');
title('����');
hold on;
grid on;
figure(5) %��ϳ̶�
plot(tt*t0,controbility(tt,1),'k-',tt*t0,controbility(tt,2),'k--','linewidth',2);
xlabel('Time (s)');
legend('���ٶȹ�ϵ','���ٶ������ؼ��ϵ');
title('��ϳ̶�');
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