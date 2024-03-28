%Çódetak
function outp=SvDelMtx(qm,k)
Tkn=T0k(qm,k);
an=Tkn(1:3,3);
pn=Tkn(1:3,4);
Van=matrix1(an);
pca=cross(pn,an);
detak=[Van,pca;0,0,0,0];
outp=detak;



