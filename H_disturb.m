function outp=H(qm)
DH=DH_Array(qm);
T01=DH2T(DH(1,:));
T12=DH2T(DH(2,:));
T23=DH2T(DH(3,:));
T34=DH2T(DH(4,:));
T45=DH2T(DH(5,:));
T56=DH2T(DH(6,:));
T67=DH2T(DH(7,:));
%基坐标系到各坐标系变换矩阵
T02=T01*T12;T03=T02*T23;
T04=T03*T34;T05=T04*T45;
T06=T05*T56;T07=T06*T67;
%基坐标系到各坐标系旋转矩阵
R01=T01(1:3,1:3);R02=T02(1:3,1:3);R03=T03(1:3,1:3);R04=T04(1:3,1:3);
R05=T05(1:3,1:3);R06=T06(1:3,1:3);R07=T07(1:3,1:3);
%各坐标系原点在基坐标系下位置
P01=T01(1:3,4);P02=T02(1:3,4);P03=T03(1:3,4);P04=T04(1:3,4);
P05=T05(1:3,4);P06=T06(1:3,4);P07=T07(1:3,4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Del1=SvDelMtx(qm,1);
Vdeta1=Del1(1:3,1:3);
deta1=Matrix1_(Vdeta1);
lamda1=Del1(1:3,4);

Del2=SvDelMtx(qm,2);
Vdeta2=Del2(1:3,1:3);
deta2=Matrix1_(Vdeta2);
lamda2=Del2(1:3,4);

Del3=SvDelMtx(qm,3);
Vdeta3=Del3(1:3,1:3);
deta3=Matrix1_(Vdeta3);
lamda3=Del3(1:3,4);

Del4=SvDelMtx(qm,4);
Vdeta4=Del4(1:3,1:3);
deta4=Matrix1_(Vdeta4);
lamda4=Del4(1:3,4);

Del5=SvDelMtx(qm,5);
Vdeta5=Del5(1:3,1:3);
deta5=Matrix1_(Vdeta5);
lamda5=Del5(1:3,4);

Del6=SvDelMtx(qm,6);
Vdeta6=Del6(1:3,1:3);
deta6=Matrix1_(Vdeta6);
lamda6=Del6(1:3,4);

Del7=SvDelMtx(qm,7);% 4*4
Vdeta7=Del7(1:3,1:3); %3*3
deta7=Matrix1_(Vdeta7);%3*1
lamda7=Del7(1:3,4);%3*1

load actual_CSof7DF
L8=zeros(3);
L7=L8+R07*Ic7*R07'+m7*MutiTwoVec(R07*Pc77+P07,R07*Pc77+P07);
L6=L7+R06*Ic6*R06'+m6*MutiTwoVec(R06*Pc66+P06,R06*Pc66+P06);
L5=L6+R05*Ic5*R05'+m5*MutiTwoVec(R05*Pc55+P05,R05*Pc55+P05);
L4=L5+R04*Ic4*R04'+m4*MutiTwoVec(R04*Pc44+P04,R04*Pc44+P04);
L3=L4+R03*Ic3*R03'+m3*MutiTwoVec(R03*Pc33+P03,R03*Pc33+P03);
L2=L3+R02*Ic2*R02'+m2*MutiTwoVec(R02*Pc22+P02,R02*Pc22+P02);
L1=L2+R01*Ic1*R01'+m1*MutiTwoVec(R01*Pc11+P01,R01*Pc11+P01);%3*3

H8=zeros(3,1);
H7=H8+m7*(R07*Pc77+P07);
H6=H7+m6*(R06*Pc66+P06);
H5=H6+m5*(R05*Pc55+P05);
H4=H5+m4*(R04*Pc44+P04);
H3=H4+m3*(R03*Pc33+P03);
H2=H3+m2*(R02*Pc22+P02);
H1=H2+m1*(R01*Pc11+P01);%3*1

U1=SvJU(L1);
U2=SvJU(L2);
U3=SvJU(L3);
U4=SvJU(L4);
U5=SvJU(L5);
U6=SvJU(L6);
U7=SvJU(L7);%3*3

X1=U1*deta1+cross(H1,lamda1);
X2=U2*deta2+cross(H2,lamda2);
X3=U3*deta3+cross(H3,lamda3);
X4=U4*deta4+cross(H4,lamda4);
X5=U5*deta5+cross(H5,lamda5);
X6=U6*deta6+cross(H6,lamda6);
X7=U7*deta7+cross(H7,lamda7);%3*1

mder8=0;
mder7=mder8+m7;
mder6=mder7+m6;
mder5=mder6+m5;
mder4=mder5+m4;
mder3=mder4+m3;
mder2=mder3+m2;
mder1=mder2+m1;

Y1=cross(deta1,H1)+mder1*lamda1;
Y2=cross(deta2,H2)+mder2*lamda2;
Y3=cross(deta3,H3)+mder3*lamda3;
Y4=cross(deta4,H4)+mder4*lamda4;
Y5=cross(deta5,H5)+mder5*lamda5;
Y6=cross(deta6,H6)+mder6*lamda6;
Y7=cross(deta7,H7)+mder7*lamda7;%3*1


D11=deta1'*X1+lamda1'*Y1;
D12=deta1'*X2+lamda1'*Y2;
D13=deta1'*X3+lamda1'*Y3;
D14=deta1'*X4+lamda1'*Y4;
D15=deta1'*X5+lamda1'*Y5;
D16=deta1'*X6+lamda1'*Y6;
D17=deta1'*X7+lamda1'*Y7;
D22=deta2'*X2+lamda2'*Y2;
D23=deta2'*X3+lamda2'*Y3;
D24=deta2'*X4+lamda2'*Y4;
D25=deta2'*X5+lamda2'*Y5;
D26=deta2'*X6+lamda2'*Y6;
D27=deta2'*X7+lamda2'*Y7;
D33=deta3'*X3+lamda3'*Y3;
D34=deta3'*X4+lamda3'*Y4;
D35=deta3'*X5+lamda3'*Y5;
D36=deta3'*X6+lamda3'*Y6;
D37=deta3'*X7+lamda3'*Y7;
D44=deta4'*X4+lamda4'*Y4;
D45=deta4'*X5+lamda4'*Y5;
D46=deta4'*X6+lamda4'*Y6;
D47=deta4'*X7+lamda4'*Y7;
D55=deta5'*X5+lamda5'*Y5;
D56=deta5'*X6+lamda5'*Y6;
D57=deta5'*X7+lamda5'*Y7;
D66=deta6'*X6+lamda6'*Y6;
D67=deta6'*X7+lamda6'*Y7;
D77=deta7'*X7+lamda7'*Y7;

H=[D11,D12,D13,D14,D15,D16,D17;D12,D22,D23,D24,D25,D26,D27;...
    D13,D23,D33,D34,D35,D36,D37;D14,D24,D34,D44,D45,D46,D47;...
    D15,D25,D35,D45,D55,D56,D57;D16,D26,D36,D46,D56,D66,D67;...
    D17,D27,D37,D47,D57,D67,D77];
outp=H;