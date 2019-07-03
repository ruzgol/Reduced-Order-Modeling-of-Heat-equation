
%% Header
% @program name: abc_example3.m
% @dependency: MATLAB symbolic toolbox required, abc_gs.m
% @task: Compute POD and GM of the convecting ABC flow
% @author: D.M. Luchtenburg
% @see: Technical Report 01/2009 by DML, BRN and MS
% @created: April 2009, DML
close all; clear all; clc;
%% Globals (defined for integration GS)
global qijk
%% Symbols
syms x y z t c1 c2 c3 A B C T w intg;
%% Parameters
% coefficients ABC flow
A = sym(sqrt(3));
B = sym(sqrt(2));
C = sym(1);
% convection velocities
c1 = sym(1);
c2 = sym(sqrt(2));
c3 = sym(sqrt(3));
T = sym(2*pi); % period for time average
nDim = 3; % number of dimensions
nModes = 6; % number of modes
eps = 1e-6; % tolerance for Runge Kutta GS integration
nu = 0; % 1 / Re
%% Init
v = zeros(nDim,nModes); v = sym(v); % temporal coefficients
b = zeros(nModes,1); b = sym(b); % spatial basis vectors
ub = zeros(nDim,1); ub = sym(ub); % base flow
u0 = zeros(nDim,1); u0 = sym(u0); % mean flow
uc = zeros(nDim,1); uc = sym(uc); % velocity field
up = zeros(nDim,1); up = sym(up); % fluctuation part velocity field
Q = zeros(nModes,nModes); Q = sym(Q); % Matrix needed for eigenprob.
%% Define base velocity field
ub = [c1; c2; c3];
%% Define basis vectors
% v_1, v_2
v(2,1) = sin(x); v(2,2) = -cos(x);
v(3,1) = cos(x); v(3,2) = sin(x);
% v_3, v_4
v(1,3) = cos(y); v(1,4) = sin(y);
v(3,3) = sin(y); v(3,4) = -cos(y);
% v_5, v_6
v(1,5) = sin(z); v(1,6) = -cos(z);
v(2,5) = cos(z); v(2,6) = sin(z);
%% Define temporal coefficients
% b_1, b_2
b(1) = B*(cos(c1*t)-sin(2*pi*c1)/(2*pi*c1));
b(2) = B*(sin(c1*t)+(cos(2*pi*c1)-1)/(2*pi*c1));
b(3) = C*(cos(c2*t)-sin(2*pi*c2)/(2*pi*c2));
b(4) = C*(sin(c2*t)+(cos(2*pi*c2)-1)/(2*pi*c2));
b(5) = A*(cos(c3*t)-sin(2*pi*c3)/(2*pi*c3));
b(6) = A*(sin(c3*t)+(cos(2*pi*c3)-1)/(2*pi*c3));
%% Construct fluctuation velocity field
for i = 1:nModes
up = up + b(i)*v(:,i);
end
%% Mean velocity field
u0 = ub + B*sin(2*pi*c1)/(2*pi*c1)*v(:,1) + B*(1-cos(2*pi*c1))/(2*pi*c1)*v(:,2) + ...
C*sin(2*pi*c2)/(2*pi*c2)*v(:,3) + C*(1-cos(2*pi*c2))/(2*pi*c2)*v(:,4) + ...
A*sin(2*pi*c3)/(2*pi*c3)*v(:,5) + A*(1-cos(2*pi*c3))/(2*pi*c3)*v(:,6);
%% Compute Q matrix
for i = 1:nModes
    for j = 1:nModes
        for k = 1:nModes
            for l = 1:nDim
                Q(i,j) = Q(i,j) + (1/T)*int(b(i)*b(k),t,0,2*pi)*...
               int( int( int( v(l,k)*v(l,j), x,0,2*pi ), y,0,2*pi ), z,0,2*pi );
            end
        end
    end
end
Q = eval(Q);
%% Compute eigenvalues
[eigenVectorArray,eigenValues] = eig(Q);
eigenValues = diag(eigenValues);
[eigenValuesTemp,I] = sort(eigenValues,'descend');
eigenValues = eigenValues(I);
eigenVectorArray = eigenVectorArray(:,I);
%% Assign weights c, u_i = sum_j c_ji v_j
c = sym(eigenVectorArray);
%% Compute POD modes
u = zeros(nDim,nModes); u = sym(u);
for i = 1:nModes
% sum
        for j = 1:nModes
            u(:,i) = u(:,i) + c(j,i)*v(:,j);
      end
end
% Compute weights and normalized modes
for i = 1:nModes
    w = int( int( int( u(:,i).'*u(:,i), x,0,2*pi ), y,0,2*pi ), z,0,2*pi );
    w = sqrt(w);
    u(:,i) = u(:,i) ./ w;
    c(:,i) = c(:,i) ./ w;
end
c = eval(c);
%% Compute Fourier coefficients
a = zeros(nModes,1); a = sym(a);
for i = 1:nModes
    a(i) = int( int( int( up.'*u(:,i), x,0,2*pi ), y,0,2*pi ), z,0,2*pi );
end
%% Compute Galerkin projection
umm = zeros(nDim,nModes+1); umm = sym(umm);
umm(:,1) = u0; umm(:,2:nModes+1) = u; % store: mean, modes
ug = zeros(nDim,nDim); ug = sym(ug); % temporary storage
qijk = zeros(nModes,nModes+1,nModes+1); qijk = sym(qijk); % convection term
lij = zeros(nModes,nModes+1); lij = sym(lij); % dissipation term
lh = ones(nDim,1); lh = sym(lh); % ones(3,1)
for i = 1:nModes
for j = 1:nModes+1 % mean flow + modes
% Dissipation term: Int dx ( (lapl u_j) * u_i )
% Compute 2nd derivatives
for ig = 1:nDim
ug(ig,1) = diff(umm(ig,j),'x',2);
ug(ig,2) = diff(umm(ig,j),'y',2);
ug(ig,3) = diff(umm(ig,j),'z',2);
end
% integrand: lapl(u_j) . u_i
intg = (ug*lh).'*umm(:,i+1); % note i-index: only POD modes!
lij(i,j) = int( int( int( intg, x,0,2*pi ), y,0,2*pi ), z,0,2*pi);
% convection term: Int dV ((grad u_k) u_j) . u_i
for k = 1:nModes+1 % mean flow + modes
% grad u_k
for ig = 1:3
ug(ig,1) = diff(umm(ig,k),'x');
ug(ig,2) = diff(umm(ig,k),'y');
ug(ig,3) = diff(umm(ig,k),'z');
end
% integrand: ((grad u_k) u_j) . u_i
intg = (ug*umm(:,j)).'*umm(:,i+1); % note i-index: only POD modes!
qijk(i,j,k) = -int( int( int( intg, x,0,2*pi ), y,0,2*pi ), z,0,2*pi );
end
end
end
% Numerical values
lij = eval(lij);
for k = 1:nModes+1 % mean flow + modes
qijk(:,:,k) = eval(qijk(:,:,k));
end
% Pressure term
% f_i = - Int grad p . u_i = Int D/Dt u . u_i
fp = zeros(nModes,1); fp = sym(fp);
dudt = zeros(nDim,1); dudt = sym(dudt);
for i = 1:nDim
dudt(i) = diff(uc(i),'t');
ug(i,1) = diff(uc(i),'x');
ug(i,2) = diff(uc(i),'y');
ug(i,3) = diff(uc(i),'z');
end
for i = 1:nModes
% integrand: D/Dt u . u_i
intg = (dudt + (ug*uc)).'*umm(:,i+1);
fp(i) = int( int( int( intg, x,0,2*pi ), y,0,2*pi ), z,0,2*pi );
end
%% Integrate Galerkin system
a0 = zeros(nModes,1);
for i = 1:nModes
a0(i) = double(subs(a(i),0)); % get initial condition
end
qijk = double(qijk);
lij = double(lij);
% Integrate with Runge Kutta
options = odeset('RelTol',eps,'AbsTol',eps*ones(nModes,1));
[tInt,aInt] = ode45(@abc_gs,[0 2*pi],a0,options);
%% Figures
close all
% Plot eigenvalues
figure
semilogy([1:nModes],(eigenValues./sum(eigenValues)),'--o')
grid on
xlabel('i'); ylabel('\lambda_i / \Sigma \lambda_i')
set(gca, 'XTick', [1:nModes]);
title('Eigenvalues')
% Plot weights
for i = 1:nModes
figure
plot([1:nModes],c(:,i),'--o')
axis([1 nModes -0.06 0.06])
xlabel('i'); ylabel('weights')
set(gca, 'XTick', [1:nModes]);
title(['Weights: T_', num2str(i),'j']);
end
% Plot Fourier coefficients
for i = 1:2:nModes
figure
p1 = ezplot(a(i),[0 eval(T)]);
set(p1,'Color','blue','LineStyle','-')
hold on
p2 = ezplot(a(i+1),[0 eval(T)]);
set(p2,'Color','red','LineStyle','--')
title('')
xlabel('t'); ylabel('a_i');
legend(['a_', num2str(i)],['a_', num2str(i+1)])
if( i==1 )
axis([0 2*pi -45 45])
elseif( i == 3 )
axis([0 2*pi -35 35])
else
axis([0 2*pi -10 10])
end
title('Fourier coefficients')
end
% Plot dissipation term lij
figure
image(lij,'CDataMapping','scaled')
axis('equal')
xlabel('j'); ylabel('i');
set(gca, 'XTickLabel', [0:nModes]);
colorbar
title('dissipation term: l_{ij}')
% Plot q_i0k (convection)
figure
image(squeeze(double(qijk(:,1,:))),'CDataMapping','scaled')
axis('equal')
cm = max(abs(min(min(squeeze(qijk(:,:,1))))),max(max(squeeze(qijk(:,:,1)))));
caxis([-cm cm])
xlabel('k'); ylabel('i');
set(gca, 'XTickLabel', [0:nModes]);
colorbar
xlabel('k'); ylabel('i');
title('q_{i0k} (->convection)')
% Plot q_ij0 (production)
figure
image(squeeze(qijk(:,:,1)),'CDataMapping','scaled')
axis('equal')
cm = max(abs(min(min(squeeze(qijk(:,:,1))))),max(max(squeeze(qijk(:,:,1)))));
caxis([-cm cm])
xlabel('j'); ylabel('i');
set(gca, 'XTickLabel', [0:nModes]);
colorbar
xlabel('j'); ylabel('i');
title('q_{ij0} (->production)')
% Plot Comparison exact Fourier coefficients and integration GS
for i = 1:2:nModes
figure
hold on
plot(tInt,aInt(:,i ),'b.-')
plot(tInt,aInt(:,i+1),'r.-')
p1 = ezplot(a(i),[0 eval(T)]);
set(p1,'Color','blue','LineStyle','-')
p2 = ezplot(a(i+1),[0 eval(T)]);
set(p2,'Color','red','LineStyle','--')
legend('a_1','a_2')
% axis
if( i==1 )
axis([0 2*pi -45 45])
elseif( i == 3 )
axis([0 2*pi -35 35])
else
axis([0 2*pi -10 10])
end
title('Comparison exact (-) Fourier coefficient and integration GS (.-)')
xlabel('t'); ylabel(['a_', num2str(i),', ','a_', num2str(i+1)]);
end
MATLAB file: abc gs.m
%% Header
% @program name: abc_gs.m
% @dependency: -
% @task: Compute GS
% @author: D.M. Luchtenburg
% @see: Technical Report 01/2009 by DML, BRN and MS
% @created: March 2009, DML
function [dadt] = abc_gs(t,a)
% compute GS: dadt = qijk*aj*ak
global qijk
nModes = length(a);
c = [1; a]; % define: a0 = 1
dadt = zeros(nModes,1);
% Sum q_ijk*a_j*a_k
for iMode = 1:nModes
for j = 1:(nModes+1)
    for k = 1:(nModes+1)
        dadt(iMode) = dadt(iMode) + qijk(iMode,j,k)*c(j)*c(k);
    end
end
end
end