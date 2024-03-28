
% theta/d/a/alpha
L(1) = Link([0 0.6 0 pi/2],'modified');
L(2) = Link([pi/2 0.5 0 -pi/2],'modified');
L(3) = Link([0 0.5 5 0],'modified');
L(4) = Link([0 0.5 5 0],'modified');
L(5) = Link([0 0.5 0 pi/2],'modified');
L(6) = Link([-pi/2 0.5 0 -pi/2],'modified');
L(7) = Link([0 0.6 0 0],'modified');

nice = SerialLink(L, 'name', 'nice');
nice.display();
nice.plot([0 pi/2 0 0 0 -pi/2 0]);%%%œ‘ æ
nice.teach;  %%% æΩÃ
