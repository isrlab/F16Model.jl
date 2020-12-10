using SymEngine

@vars x y z phi theta psi vt alpha beta P Q R

st    = sin(theta);
ct    = cos(theta);
tt    = tan(theta);

sphi  = sin(phi);
cphi  = cos(phi);
spsi  = sin(psi);
cpsi  = cos(psi);

sa    = sin(alpha); # sin(alpha) 
ca    = cos(alpha); # cos(alpha) 
sb    = sin(beta); # sin(beta)  
cb    = cos(beta); # cos(beta)  
tb    = tan(beta); # tan(beta)  


# Transformation Matrix
T1 = [1      0           0; 
      0  cos(phi) sin(phi); 
      0 -sin(phi) cos(phi)];

T2 = [cos(theta) 0 -sin(theta); 
       0         1        0; 
      sin(theta) 0 cos(theta)];

T3 = [cos(psi) sin(psi) 0;
     -sin(psi) cos(psi) 0;
         0       0      1];

Tbo = T1*T2*T3;
Tob = transpose(Tbo);

U = vt*ca*cb;  # directional velocities.
V = vt*sb;
W = vt*sa*cb;

xd = U*(ct*cpsi) + V*(sphi*cpsi*st - cphi*spsi) + W*(cphi*st*cpsi + sphi*spsi); # npos dot
yd = U*(ct*spsi) + V*(sphi*spsi*st + cphi*cpsi) + W*(cphi*st*spsi - sphi*cpsi); # epos dot
zd = U*st - V*(sphi*ct) - W*(cphi*ct); # alt dot

Xd = Tob*[U;V;W];

## Kinematic equations
phid = P + tt*(Q*sphi + R*cphi); # phi dot
thd = Q*cphi - R*sphi; # theta dot
psid = (Q*sphi + R*cphi)/ct; # psi dot

f = [xd,yd,zd,phid,thd,psid];
X = [x,y,z,phi,theta,psi,vt,alpha,beta,P,Q,R];
X0 = [0,0,0,0,0,0,300,0,0,0,0,0];

function jacobian(f,X,X0)::Tuple{Matrix,Matrix}
    n = length(f);
    m = length(X);
    A = zeros(Basic,n,m);
    As = zeros(Basic,n,m);
    for i=1:n
        for j=1:m
            a = diff(f[i],X[j]);
            As[i,j] = a;
            for k=1:m
                a = subs(a,X[k],X0[k]);
            end
            A[i,j] = a;
        end
    end
    return (As,A);
end

As, A = jacobian(f,X,X0);