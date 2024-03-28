%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%     ���ƹ滮     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% ����ѧϰ���ƹ滮  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% ��������Ŀ��ƹ滮  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%     �����Ŷ�     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc;
load actual_CSof7DF;
tic
t0=0.05; %ʱ�䲽��
n1=150; %�����ؽڹ滮����
n=700; %�����ؽڹ滮�������ܲ���
N=100; %����ѭ������
wn=1;
sigma=0.9;
kp=wn^2;
kd=2*sigma*wn;
miu=9.8; %0,1֮�䣿�����ٶ����֮�����ü��ٶȽ�����ܴ����ٶ�

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%    ��ʼ����    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% qint=[-10,120,150,30,20,-120,50]'*pi/180;%��ʼ�ؽڽ�
% qint=[50,150,60,-60,130,-170,30]'*pi/180;%��ʼ�ؽڽ�
qint=[40,150,60,-60,15,-85,40]'*pi/180;%��ʼ�ؽڽ�
q=zeros(7,n+1);
q(:,1)=qint;
vq=zeros(7,n+1);
aq=zeros(7,n);
taom=zeros(7,n);
failnum=4;  %���Ϲؽ�
activenum=3; %�����ؽ�
q_pd=47*pi/180;vq_pd=0;  %�����ؽ������Ƕȡ����ٶ�
qd=zeros(7,1);
qd(:,1)=qint;
delta=15*pi/180;qd(activenum,1)=qint(activenum,1)+delta; %��һ��ѭ����ʼ�����ؽ������Ƕ�
error=zeros(N,2);%����ÿ�ε���ѭ�������������
for j=1:N
j
for i=1:n1   %%�ؽڵĹ켣�滮
    aq(:,i)=kp*(qd(:,1)-q(:,i))-kd*vq(:,i);   %PD������
    vq(:,i+1)=vq(:,i)+aq(:,i)*t0;
    q(:,i+1)=q(:,i)+vq(:,i)*t0;
end
%%%�����ؽڵĹ켣�滮
a0=qint(activenum,1);a1=0;a2=0;a3=10*(qd(activenum,1)-qint(activenum,1))/(t0*n1)^3;
a4=-15*(qd(activenum,1)-qint(activenum,1))/(t0*n1)^4;a5=6*(qd(activenum,1)-qint(activenum,1))/(t0*n1)^5;
for i=1:n1
    aq(activenum,i)=2*a2+6*a3*(t0*i)+12*a4*(t0*i)^2+20*a5*(t0*i)^3;
    vq(activenum,i)=a1+2*a2*(t0*i)+3*a3*(t0*i)^2+4*a4*(t0*i)^3+5*a5*(t0*i)^4;
    q(activenum,i)=a0+a1*(t0*i)+a2*(t0*i)^2+a3*(t0*i)^3+a4*(t0*i)^4+a5*(t0*i)^5;
end
for i=n1+1:n    %%%�ؽ�ֹͣ���������
    q(:,i)=q(:,n1);
end

%%%�����ؽڵĹ켣�滮
for i=1:n   
    D=H_disturb(q(:,i));
    C=FCK_disturb(q(:,i),vq(:,i));
    aq(failnum,i)=(-miu*vq(failnum,i)-C(failnum,1)-(D(failnum,1)*aq(1,i)+D(failnum,2)*aq(2,i)+D(failnum,3)*aq(3,i)+...
        D(failnum,5)*aq(5,i)+D(failnum,6)*aq(6,i)+D(failnum,7)*aq(7,i)))/D(failnum,failnum);
    vq(failnum,i+1)=vq(failnum,i)+aq(failnum,i)*t0;
    q(failnum,i+1)=q(failnum,i)+vq(failnum,i)*t0;
    taom(:,i)=H(q(:,i))*aq(:,i)+FCK(q(:,i),vq(:,i));
    
    e1=q_pd-q(failnum,i+1);
    e2=vq_pd-vq(failnum,i+1);
    if (i>n1) && (vq(failnum,i+1) < 1.5e-1*pi/180) %%%����ʱ�ȸĵ���ֵ��������
        break;
    end
end
if i<=n-1
    for k=i+1:n+1
        q(failnum,k+1)=q(failnum,i+1);
        vq(failnum,k+1)=vq(failnum,i+1);
    end
end

%%%�жϵ���
e1*180/pi
e2*180/pi
error(j,1)=e1; %�Ƕ����
error(j,2)=e2; %�ٶ����
if (abs(e1)<=1e-1*pi/180)   %�����������
    break;
else if  error(j,1) < (q_pd-qd(failnum,1))           %%%�ص���Ʒ�����
        tp=0.4; %%k����0,5������]����ֵ�Ե���������Ӱ��
%         k=delta*(1+error(j,1)/(q(failnum,i+1)-q(failnum,1)));
%         qd(activenum,1) = qd(activenum,1)+k1*k+k2*error(j,2);   %�������������������ٶȵķ�����Ӧ�÷����������ؽڵ������ٶȣ����������нǶȣ���
%         if j==1
            qd(activenum,1) = qd(activenum,1)+tp*error(j,1)*sign(q(activenum,n1)-q(activenum,1));
%             if (abs(e1)<=1e-3) && (abs(e2)>1e-4)
%                 n1=n1+10;
%             end
%         else
%             qd(activenum,1) = qd(activenum,1)+k1*error(j,1)*sign(q(activenum,n1)-q(activenum,1))+k2*(error(j,1)-error(j-1,1));
%         end
    else
        qd(activenum,1) = qd(activenum,1)-2*delta;  %%�����ؽڷ����˶�
    end
end
end

%%%%% �ؽ�����(��) %%%%%%
tt=1:n;
figure(1) %����ѭ�������������
s=1:20;
plot(s,abs(error(s,1)*180/pi),'k-','linewidth',2);
xlabel('k');
ylabel('Output Error (��)');
axis([1 10 0 110]);
set(gca,'xtick',1:10);
grid on;
hold on;
%%%�����Ƕ�
desire_qp=zeros(n,1);
for i=1:n
    desire_qp(i,1)=47;
end
figure(2) %�Ƕ�
plot(tt*t0,q(activenum,tt)*180/pi,'k-',tt*t0,q(failnum,tt)*180/pi,'r-',tt*t0,desire_qp(tt,1),'r--','linewidth',2);
xlabel('Time (s)');
ylabel('Angle (��)');
legend('The Active Joint','The Passive Joint','The Desire Angle');
hold on;
grid on;
figure(3) %���ٶ�
plot(tt*t0,vq(activenum,tt)*180/pi,'k-',tt*t0,vq(failnum,tt)*180/pi,'r-','linewidth',2);
xlabel('Time (s)');
ylabel('Vlocity (��/s)');
legend('The Active Joint','The Passive Joint');
hold on;
grid on;
figure(4) %�Ǽ��ٶ�
plot(tt*t0,aq(activenum,tt)*180/pi,'k-',tt*t0,aq(failnum,tt)*180/pi,'r-','linewidth',2);
xlabel('Time (s)');
ylabel('Acceleration (��/s^2)');
legend('The Active Joint','The Passive Joint');
hold on;
grid on;

toc