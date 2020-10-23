clear all;
syms ds dphi dds ddphi s phi l m1 m2 I1 I2 g real; 

Zzero = [0 0 0]';
Z0 = [0 0 1]';
Z1 = [1 0 0]';

T_1 = simplify(Rz(phi)*Tx(s+l));

px1 = T_1(1,4);
py1 = T_1(2,4);
pz1 = T_1(3,4);

Jv1 = [diff(px1, phi), diff(px1, s);
      diff(py1, phi), diff(py1, s);
      diff(pz1, phi), diff(pz1, s)];
  
Jv2 = [Z1, Z0];

Jw1 = [Z0, Zzero];
Jw2 = [0];

Rz = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];

D1 = m1*Jv1'*Jv1 + Jw1'*Rz*(I1)*Rz'*Jw1;

D2 = m2*Jv2'*Jv2 + 0; 

D = D1+D2;
D = simplify(D)

P1 = m1*g*0;
P2 = m2*g*((l+s)*sin(phi));
P = P1+P2
G1 = diff(P, phi);
G2 = diff(P, s);

G = [G1; G2]

q = [phi; s];
dq = [dphi; ds];
ddq = [ddphi; dds];
C = Coriolis(D,q,dq,2);
C=simplify(C)

tor = D*ddq+C*dq+G
D(phi,s) = subs(D,{m1, m2, I1, I2, l},{2 2 2 1 0.2});
C(phi,s,dphi,ds) = subs(C*dq ,{m1, m2, I1, I2, l},{2 2 2 1 0.2});
G(phi,s) = subs(G ,{m1, m2, I1, I2, l, g},{2 2 2 1 0.2 9.81});

q1_0 = 0;
q2_0 = 0;
dq1_0 = 0;
dq2_0 = 0;
dt = 0.01;
U = [0;0];
n=1000;
for i = 1:n

 q1p(i) = q1_0;
 q2p(i) = q2_0;
 dq1p(i) = dq1_0;
 dq2p(i) = dq2_0;
 ddq = inv(D(q1_0, q2_0))*(U-C(q1_0, q2_0,dq1_0,dq2_0)-G(q1_0,q2_0));

 dq1_0 = dq1p(i) + double(ddq(1)*dt);
 dq2_0 = dq2p(i) + double(ddq(2)*dt);

 q1_0 = q1p(i) + dq1_0*dt;
 q2_0 = q2p(i) + dq2_0*dt;
end
t = 0:0.1:(0.1*(n-1));
figure
plot(t,dq1p)

figure
plot(t,dq2p)


