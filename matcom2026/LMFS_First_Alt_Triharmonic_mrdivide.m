% First alternative approach    
% \Delta^3 u(x,y,z)=0, (x,y,z)\in Omega
% BC: Dirichlet, Neumann, Laplacian
clear; 
u=@(x,y,z) (sqrt((x-1).^2+(y-2).^2+(z-3).^2)).^3 ;
ux=@(x,y,z) 3*(x-1).*(sqrt((x-1).^2+(y-2).^2+(z-3).^2));
uy=@(x,y,z) 3*(y-2).*(sqrt((x-1).^2+(y-2).^2+(z-3).^2));
uz=@(x,y,z) 3*(z-3).*(sqrt((x-1).^2+(y-2).^2+(z-3).^2));
uL=@(x,y,z) 12*(sqrt((x-1).^2+(y-2).^2+(z-3).^2));
DMatrix=@(A,B) sqrt( (A(:,1)-B(:,1)').^2+(A(:,2)-B(:,2)').^2+(A(:,3)-B(:,3)').^2);
R=3;      % Radius of source sphere
k=1;      % LNSS : parameter for matrix sparsification
ni=10000; % # of interior collocation points 
rh=0.2;   % radius of hollow sphere
NBH=200;  % # boundary points on every sphere
NBC=20;   % # of boundary points on every edge of cube
ns=50;    % # of collocation points in each local domain
nsh=45;   % # of source points in each local domain
% Generation of interior, boundary points, normal vectors, and source points
[xi,yi,zi,bx,by,bz,xn,yn,zn,xs,ys,zs]=cube_hollow(ni,NBC,NBH,rh,nsh);
p=haltonset(1); q=net(p,nsh);  
noise=0.3; 
% Source points on a unit sphere with noise disturbance
xs=xs.*(1+q*noise); ys=ys.*(1+q*noise); zs=zs.*(1+q*noise); 
x=[bx;xi]; y=[by;yi]; z=[bz;zi];  % All the collocation points: boundary + interior points
nb=length(bx); n=nb+ni;
idx=knnsearch([x,y,z],[x,y,z],'k',k*ns+1);  
w = zeros(n,ns); wx=zeros(nb,ns); wy=wx;wz=wx; wL=wx;
% Building local matrix at each local domain
for i=1:n
    % [tpx,tpy,tpz]: Collocation points in the local domain    
    tpx=x(idx(i,2:k:k*ns+1)); tpy=y(idx(i,2:k:k*ns+1)); tpz=z(idx(i,2:k:k*ns+1));
    sx=xs*R+x(i); sy=ys*R+y(i); sz=zs*R+z(i); % Source points
    DM=DMatrix([tpx,tpy,tpz],[sx,sy,sz]);
    G=DM.^3;
    DM1=DMatrix([x(i),y(i),z(i)],[sx,sy,sz]);
    if i>nb % Interior points
       mfs=DM1.^3; 
       w(i,:)=mfs/G; % Dirichlet (1st) BC
    else       
        tax=3*(x(i)-sx').*DM1; 
        tay=3*(y(i)-sy').*DM1;     
        taz=3*(z(i)-sz').*DM1;          
        wx(i,:)=tax/G; % Neumann (2nd) BC
        wy(i,:)=tay/G; 
        wz(i,:)=taz/G; 
        mfsL=12*DM1;  
        wL(i,:)=mfsL/G; % Laplacian (3rd) BC 
    end
end

% Establish global sparse matrix
C1=zeros(2*nb*ns+ni*(ns+1)+nb,1); C2=C1; C3=C1;
C1(1:nb)=(1:nb); C2(1:nb)=(1:nb); C3(1:nb)=1; % Dirichlet BC
ic=nb;
for i=1:nb % Neumann BC
    C1(ic+1:ic+ns)=i+nb; 
    C2(ic+1:ic+ns)=idx(i,2:k:k*ns+1);
    C3(ic+1:ic+ns)=xn(i,1)*wx(i,:)+yn(i,1)*wy(i,:)+ zn(i,1)*wz(i,:);
    ic=ic+ns;
end
for i=1:nb % Laplacian BC
    C1(ic+1:ic+ns)=i+2*nb; 
    C2(ic+1:ic+ns)=idx(i,2:k:k*ns+1);
    C3(ic+1:ic+ns)=wL(i,:);
    ic=ic+ns;
end
for i=nb+1:n % Interior 
    C1(ic+1:ic+ns+1)=i+2*nb;  
    C2(ic+ns+1)=i; C3(ic+ns+1)=1;
    C2(ic+1:ic+ns)=idx(i,2:k:k*ns+1);
    C3(ic+1:ic+ns)=-w(i,:);
    ic=ic+ns+1;
end

A=sparse(C1,C2,C3,n+2*nb,n); % Global sparse matrix
v=u(x,y,z); % Exact solution
vx=ux(bx,by,bz); vy=uy(bx,by,bz); vz=uz(bx,by,bz); 
b=zeros(n+2*nb,1);
b(1:nb)=v(1:nb); % Dirchlet BC
b(nb+1:2*nb)=xn.*vx+yn.*vy+zn.*vz; % Neumann BC
b(2*nb+1:3*nb)=uL(bx,by,bz); % Laplacian BC
app=A\b; % Approximate solution
ERROR=max(abs(app-v))/max(abs(v)); % Relative max. error
fprintf('ni = %5d, nb = %3d, ns = %2d, R = %3.1f, k = %1d \n',ni,nb,ns,R,k); 
fprintf('Error = %8.3e\n',ERROR); 