%反对称矩阵对应的向量
%与函数matrix1.m的求解过程正好相反
function outp=Matrix1_(Mt)
outp=[(Mt(3,2)-Mt(2,3))/2,(Mt(1,3)-Mt(3,1))/2,(Mt(2,1)-Mt(1,2))/2]';