function [B,C,P,T,U,R,R2X,R2Y]=simpls(X,Y,A,S,XtX);
%	      Implementação completa da abordagem SIMPLS para regressão PLS
% para (multivariada) Y.			    
%									    
% Reference.: S.de Jong, Chemom.Intell.Lab.Syst.,18 (1993) 251-263.   
% Entrada    X   centrado na média 
%            Y    centrado na média  
%            A número máximo de Variaveis latentes    
%            S = X'*Y  
%	      XtX  = X'*X 



[n,px] = size(X); [n,m] = size(Y);   			
if nargin<5, S = []; end, if isempty(S), S=(Y'*X)'; end	
if nargin<4, XtX=[]; end		
if isempty(XtX) & n>3*px, XtX = X'*X; end		
if nargin<3, A=10; end, A = min([A px n-1]);
T = zeros(n ,A); U = T;				
R = zeros(px,A); P = R; V = R;
C = zeros(m ,A); 
R2Y = zeros(1,A);
z = zeros(m,1); v = zeros(px,1);
if n>px, S0 = S; end, StS = S'*S;
nm1 = n-1;
tol = 0;
for a = 1:A
  StS = StS-z*z'; 
  [Q,LAMBDA] = eig(StS); 
  [lambda,j] = max(diag(LAMBDA)); 
  q = Q(:,j(1));
  r = S*q;
  t = X*r;
  if isempty(XtX), p = (t'*X)'; else p = XtX*r; end
  if n>px, d = sqrt(r'*p/nm1); else d = sqrt(t'*t/nm1); end
  if d<tol, 
	disp(' ')
        disp('AVISO: o número necessário de fatores (A) é muito alto!')
	disp('Menos fatores de PLS foram extraídos dos dados do programa PLSSIM!') 
	disp(' ')
	break,
, 	else tol=max(tol,d/1e5);
  end
  v = p-V(:,1:max(1,a-1))*(p'*V(:,1:max(1,a-1)))'; v = v/sqrt(v'*v); 
  z = (v'*S)'; 
  S = S-v*z'; 
							
  V(:,a) = v;
  R(:,a) = r/d; 						
  P(:,a) = p/(d*nm1); 					
  T(:,a) = t/d;							
  U(:,a) = Y*q;						
  C(:,a) = q*(lambda(1)/(nm1*d)); 			
  R2Y(1,a) =  lambda(1)/d;			
end
clear StS V LAMBDA Q p q r t v z;
if d<tol,
 A=a-1; a=A; T=T(:,1:A); U=U(:,1:A); R=R(:,1:A); P=P(:,1:A); C=C(:,1:A);
end
while a>1
  U(:,a) = U(:,a)-T(:,1:a-1)*(U(:,a)'*T(:,1:a-1)/nm1)'; 
  a=a-1; 
end
B = R*C';							% Coeficientes de regressão
if isempty(XtX), sumX2=sum(X.^2); else sumX2 = sum(diag(XtX)); end
R2X = 100*nm1/sum(sumX2)*cumsum(sum(P.^2)); 
R2Y = 100/nm1/sum(sum(Y.^2))*cumsum(R2Y(1:A).^2);



