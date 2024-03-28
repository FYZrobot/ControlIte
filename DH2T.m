%由DH参数求相邻连杆坐标系的变换矩阵
%输入依次为两z轴夹角，两z轴距离，两x轴夹角，两x轴距离,即alpha,a,theta,d
function zans=DH2T(Ar)
T=[cos(Ar(3)),-sin(Ar(3)),0,Ar(2);...
    sin(Ar(3))*cos(Ar(1)),cos(Ar(3))*cos(Ar(1)),-sin(Ar(1)),-Ar(4)*sin(Ar(1));...
    sin(Ar(3))*sin(Ar(1)),cos(Ar(3))*sin(Ar(1)),cos(Ar(1)),Ar(4)*cos(Ar(1));0,0,0,1];
zans=T;