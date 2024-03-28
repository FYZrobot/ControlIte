function outp=Svqq(a,k)
if k==1
    J=[a(1)*a(2),a(1)*a(3),a(1)*a(4),a(1)*a(5),a(1)*a(6),a(1)*a(7)]';
end
if k==2
    J=[a(2)*a(3),a(2)*a(4),a(2)*a(5),a(2)*a(6),a(2)*a(7)]';
end
if k==3
    J=[a(3)*a(4),a(3)*a(5),a(3)*a(6),a(3)*a(7)]';
end
if k==4
    J=[a(4)*a(5),a(4)*a(6),a(4)*a(7)]';
end
if k==5
    J=[a(5)*a(6),a(5)*a(7)]';
end
if k==6
    J=a(6)*a(7);
end
outp=J;