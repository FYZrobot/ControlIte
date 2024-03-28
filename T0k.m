%惯性系到{k}坐标系的变换矩阵，为求detak做准备（用于函数SvDelMtx）
function outp=T0k(qm,k)
DH=DH_Array(qm);
T01=DH2T(DH(1,:));
T12=DH2T(DH(2,:));
T23=DH2T(DH(3,:));
T34=DH2T(DH(4,:));
T45=DH2T(DH(5,:));
T56=DH2T(DH(6,:));
T67=DH2T(DH(7,:));
if k==1
    T0k=T01;
end
if k==2
    T0k=T01*T12;
end
if k==3
    T0k=T01*T12*T23;
end
if k==4
    T0k=T01*T12*T23*T34;
end
if k==5
    T0k=T01*T12*T23*T34*T45;
end
if k==6
    T0k=T01*T12*T23*T34*T45*T56;
end
if k==7
    T0k=T01*T12*T23*T34*T45*T56*T67;
end
outp=T0k;