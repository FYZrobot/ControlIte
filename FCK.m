%离心力和科氏力
function outp=FCK(qm,vqm)
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

% if ww>0.5
    load CSof7DF
% else
%     load actual_CSof7DF
% end

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

V1=L1*Vdeta1'+MutiTwoVec(H1,lamda1);
V2=L2*Vdeta2'+MutiTwoVec(H2,lamda2);
V3=L3*Vdeta3'+MutiTwoVec(H3,lamda3);
V4=L4*Vdeta4'+MutiTwoVec(H4,lamda4);
V5=L5*Vdeta5'+MutiTwoVec(H5,lamda5);
V6=L6*Vdeta6'+MutiTwoVec(H6,lamda6);
V7=L7*Vdeta7'+MutiTwoVec(H7,lamda7);%3*3

W1=SvJW(V1);
W2=SvJW(V2);
W3=SvJW(V3);
W4=SvJW(V4);
W5=SvJW(V5);
W6=SvJW(V6);
W7=SvJW(V7);%3*3

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

Z12=W2*deta1+cross(Y2,lamda1);
Z13=W3*deta1+cross(Y3,lamda1);
Z14=W4*deta1+cross(Y4,lamda1);
Z15=W5*deta1+cross(Y5,lamda1);
Z16=W6*deta1+cross(Y6,lamda1);
Z17=W7*deta1+cross(Y7,lamda1);
Z23=W3*deta2+cross(Y3,lamda2);
Z24=W4*deta2+cross(Y4,lamda2);
Z25=W5*deta2+cross(Y5,lamda2);
Z26=W6*deta2+cross(Y6,lamda2);
Z27=W7*deta2+cross(Y7,lamda2);
Z34=W4*deta3+cross(Y4,lamda3);
Z35=W5*deta3+cross(Y5,lamda3);
Z36=W6*deta3+cross(Y6,lamda3);
Z37=W7*deta3+cross(Y7,lamda3);
Z45=W5*deta4+cross(Y5,lamda4);
Z46=W6*deta4+cross(Y6,lamda4);
Z47=W7*deta4+cross(Y7,lamda4);
Z56=W6*deta5+cross(Y6,lamda5);
Z57=W7*deta5+cross(Y7,lamda5);
Z67=W7*deta6+cross(Y7,lamda6);%3*1

D112=deta1'*Z12;% 
D113=deta1'*Z13;
D114=deta1'*Z14;
D115=deta1'*Z15;
D116=deta1'*Z16;
D117=deta1'*Z17;
D122=deta2'*Z12;
D123=deta2'*Z13;
D124=deta2'*Z14;
D125=deta2'*Z15;
D126=deta2'*Z16;
D127=deta2'*Z17;
D133=deta3'*Z13;
D134=deta3'*Z14;
D135=deta3'*Z15;
D136=deta3'*Z16;
D137=deta3'*Z17;
D144=deta4'*Z14;
D145=deta4'*Z15;
D146=deta4'*Z16;
D147=deta4'*Z17;
D155=deta5'*Z15;
D156=deta5'*Z16;
D157=deta5'*Z17;
D166=deta6'*Z16;
D167=deta6'*Z17;
D177=deta7'*Z17;

D213=deta1'*Z23;
D214=deta1'*Z24;
D215=deta1'*Z25;
D216=deta1'*Z26;
D217=deta1'*Z27;
D223=deta2'*Z23;
D224=deta2'*Z24;
D225=deta2'*Z25;
D226=deta2'*Z26;
D227=deta2'*Z27;
D233=deta3'*Z23;
D234=deta3'*Z24;
D235=deta3'*Z25;
D236=deta3'*Z26;
D237=deta3'*Z27;
D244=deta4'*Z24;
D245=deta4'*Z25;
D246=deta4'*Z26;
D247=deta4'*Z27;
D255=deta5'*Z25;
D256=deta5'*Z26;
D257=deta5'*Z27;
D266=deta6'*Z26;
D267=deta6'*Z27;
D277=deta7'*Z27;

D314=deta1'*Z34;
D315=deta1'*Z35;
D316=deta1'*Z36;
D317=deta1'*Z37;
D324=deta2'*Z34;
D325=deta2'*Z35;
D326=deta2'*Z36;
D327=deta2'*Z37;
D334=deta3'*Z34;
D335=deta3'*Z35;
D336=deta3'*Z36;
D337=deta3'*Z37;
D344=deta4'*Z34;
D345=deta4'*Z35;
D346=deta4'*Z36;
D347=deta4'*Z37;
D355=deta5'*Z35;
D356=deta5'*Z36;
D357=deta5'*Z37;
D366=deta6'*Z36;
D367=deta6'*Z37;
D377=deta7'*Z37;

D415=deta1'*Z45;
D416=deta1'*Z46;
D417=deta1'*Z47;
D425=deta2'*Z45;
D426=deta2'*Z46;
D427=deta2'*Z47;
D435=deta3'*Z45;
D436=deta3'*Z46;
D437=deta3'*Z47;
D445=deta4'*Z45;
D446=deta4'*Z46;
D447=deta4'*Z47;
D455=deta5'*Z45;
D456=deta5'*Z46;
D457=deta5'*Z47;
D466=deta6'*Z46;
D467=deta6'*Z47;
D477=deta7'*Z47;% D477=deta6*Z47;

D516=deta1'*Z56;
D517=deta1'*Z57;
D526=deta2'*Z56;
D527=deta2'*Z57;
D536=deta3'*Z56;
D537=deta3'*Z57;
D546=deta4'*Z56;
D547=deta4'*Z57;
D556=deta5'*Z56;
D557=deta5'*Z57;
D566=deta6'*Z56;
D567=deta6'*Z57;
D577=deta7'*Z57;

D617=deta1'*Z67;
D627=deta2'*Z67;
D637=deta3'*Z67;
D647=deta4'*Z67;
D657=deta5'*Z67;
D667=deta6'*Z67;
D677=deta7'*Z67;

MB1=[D112, D113, D114, D115, D116, D117; 0, D213, D214, D215, D216, D217;...
    -D213, 0, D314, D315, D316, D317; -D214, -D314, 0, D415, D416, D417;...
    -D215, -D315, -D415, 0, D516, D517; -D216, -D316, -D416, -D516, 0, D617;...
    -D217, -D317, -D417, -D517, -D617, 0];%7*6
MB2=[D123, D124, D125, D126, D127;D223, D224, D225, D226, D227;...
    0,D324, D325, D326, D327;-D324, 0, D425, D426, D427;...
    -D325, -D425, 0, D526, D527;-D326, -D426, -D526, 0, D627;...
    -D327, -D427, -D527,-D627, 0];%7*5
MB3=[D134, D135, D136, D137; D234, D235, D236, D237; D334, D335,D336, D337;...
    0, D435, D436, D437; -D435, 0, D536, D537; -D436,-D536, 0, D637;...
    -D437, -D537, -D637, 0];%7*4
MB4=[D145, D146, D147; D245, D246, D247; D345, D346, D347; D445, D446, D447;...
    0, D546, D547; -D546, 0, D647; -D547, -D647, 0];%7*3
MB5=[D156, D157; D256, D257; D356, D357; D456, D457; D556, D557; 0, D657; -D657, 0];%7*2
MB6=[D167; D267; D367; D467; D567; D667; 0];%7*1
MC=[0, D122, D133, D144, D155, D166, D177; -D112, 0, D233, D244, D255, D266, D277; ...
    -D113, -D223, 0, D344, D355, D366, D377; -D114, -D224, -D334, 0, D455, D466, D477;...
    -D115, -D225, -D335, -D445, 0, D566, D577; -D116, -D226, -D336, -D446, -D556, 0, D677;...
    -D117, -D227, -D337, -D447, -D557, -D667, 0];%7*7
FK=2*(MB1*Svqq(vqm,1)+MB2*Svqq(vqm,2)+MB3*Svqq(vqm,3)...
    +MB4*Svqq(vqm,4)+MB5*Svqq(vqm,5)+MB6*Svqq(vqm,6));%科氏力项
FC=MC*SvVec2(vqm);%离心力项
F=FK+FC;
outp=F;

