%��DH������������������ϵ�ı任����
%��������Ϊ��z��нǣ���z����룬��x��нǣ���x�����,��alpha,a,theta,d
function zans=DH2T(Ar)
T=[cos(Ar(3)),-sin(Ar(3)),0,Ar(2);...
    sin(Ar(3))*cos(Ar(1)),cos(Ar(3))*cos(Ar(1)),-sin(Ar(1)),-Ar(4)*sin(Ar(1));...
    sin(Ar(3))*sin(Ar(1)),cos(Ar(3))*sin(Ar(1)),cos(Ar(1)),Ar(4)*cos(Ar(1));0,0,0,1];
zans=T;