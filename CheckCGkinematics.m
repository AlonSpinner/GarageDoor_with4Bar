Rx=real(rc);
Ry=imag(rc);

Vx=real(DoorCGvel);
Vy=imag(DoorCGvel);

Ax=real(DoorCGacce);
Ay=imag(DoorCGacce);

VxNum=diff(Rx)./diff(Tvec);
VyNum=diff(Ry)./diff(Tvec);

AxNum=diff(VxNum)./diff(Tvec(1:end-1));
AyNum=diff(VyNum)./diff(Tvec(1:end-1));

plot(Tvec(1:end-2),AyNum);
hold on
plot(Tvec,Ay);


