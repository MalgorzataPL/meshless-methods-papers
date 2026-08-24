%Points generator in a cube with holes
%NI: # of interior points
%NBC: # of boundary points on every edge of cube 
%NBH: # boundary points on every sphere
%rh:    radius of each hollow sphere
%nsh: # of source points in each local domain
%[xi,yi,zi]: interior points
%[xb,yb,zb]: boundary points
%[xn,yn,zn]: normal vectors on the boundary
%[xs,ys,zs]: source points on a unit sphere
function [xi,yi,zi,xb,yb,zb,xn,yn,zn,xs,ys,zs]=cube_hollow(NI,NBC,NBH,rh,nsh)
bdpt=sphere_ran(NBH);  %reference sphere
normalvec0=bdpt;
x0=bdpt(:,1); y0=bdpt(:,2);  z0=bdpt(:,3);
%%%%%% setting down the sphere centres in cube
gap=0.15; % gap from cube edges
n=(2-2*gap)/(2*rh); % max # of sphere on every direction
a = linspace(-1-rh,1+rh,floor(n)+1);
[X,Y,Z] = meshgrid(a(2:end-1)); %sphere of hollow
xc=X(:); yc=Y(:);zc=Z(:); % sphere centers
NH=length(xc); 
xbh=zeros(NBH,NH); ybh=xbh; zbh=xbh;
for i=1:length(xc)
     xbh(:,i)=rh*x0+xc(i);
     ybh(:,i)=rh*y0+yc(i);
     zbh(:,i)=rh*z0+zc(i);
end
xbh=xbh(:); ybh=ybh(:); zbh=zbh(:);
normalvec=-repmat(normalvec0,NH,1);% direction should go inside
xbhn=normalvec(:,1);ybhn=normalvec(:,2);zbhn=normalvec(:,3);
%% interior points picking
pp=haltonset(3,'skip',200); p=net(pp,5*NI);
xp=2*p(:,1)-1;  yp=2*p(:,2)-1;  zp=2*p(:,3)-1;
j=1; xi1=zeros(5*NI,1); yi1=xi1;zi1=xi1;
for ii=1:length(xp)
     d=sqrt((xp(ii)-xc).^2+(yp(ii)-yc).^2+(zp(ii)-zc).^2);
      if d>rh
          xi1(j)=xp(ii);  yi1(j)=yp(ii);   zi1(j)=zp(ii);
           j=j+1;
      end
end
xi=xi1(1:NI);yi=yi1(1:NI);zi=zi1(1:NI);
%% cube boundary
x1=-1;x2=1;y1=-1;y2=1;z1=-1;z2=1; 
a1=linspace(-1,1,NBC);  
[X,Y,Z]=ndgrid(a1); 
Px=X(:); Py=Y(:); Pz=Z(:);
id = Px==x1 |Px==x2 |Py==y1 |Py==y2 |Pz==z1 |Pz==z2;
xbc=Px(id);  ybc=Py(id);   zbc=Pz(id);
%% cube normal
p=[xbc, ybc, zbc];
a = 1; % size of cubic
n=[0 0 0];
for iii=1:size(p,1)
    for jj=1:3
        if (p(iii,jj)==a)
             n(iii,jj)=1;
         elseif (p(iii,jj)==-a)
             n(iii,jj)=-1;
        end    
    end
n(iii,:)=n(iii,:)./sqrt(n(iii,1)^2+n(iii,2)^2+n(iii,3)^2);
end
xbcn=n(:,1);ybcn=n(:,2);zbcn=n(:,3);
%% collecting boundary & its normal
xb=[xbc;xbh];   yb=[ybc;ybh];   zb=[zbc;zbh];
xn=[xbcn;xbhn]; yn=[ybcn;ybhn]; zn=[zbcn;zbhn];
%% Source points on a unit sphere
sources=sphere_ran(nsh);
xs=sources(:,1); ys=sources(:,2); zs=sources(:,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%
function bdpt=sphere_ran(n)
bdpt=zeros(n,3);
dlong = pi*(3-sqrt(5));
long = 0;
dz = sqrt(1)*2.0/n;
z = sqrt(1) - dz/2;
for k = 0 : n-1
   r = sqrt(1-z.^2);
   bdpt(k+1,2)=cos(long)*r;
   bdpt(k+1,3)=sin(long)*r;
   bdpt(k+1,1)=z;
   z = z - dz;
   long = long + dlong;
end