function [model] = my_calc_qt_limits_val (Xcal,Ycal,Xnew,nvl);

	[Xce,mX,Xnew_ce]=centrar(Xcal,1,Xnew);		
	[Yce,mY]=centrar(Ycal,1);
	[B,C,P,T,U,R,R2X,R2Y]=simpls(Xce,Yce,nvl,Xce'*Yce,[]);	        
	Yhat=Xce*B;
    
	[Yhat]=centrar_inverso(Yhat,mY);
    
    
    Ynew_hat=Xnew_ce*B;			
	[Ynew_hat]=centrar_inverso(Ynew_hat,mY);
    
    
    Xrecuperado = T*P';
    [Xrecuperado]=centrar_inverso(Xrecuperado,mX);
    Qcont=Xcal-Xrecuperado;
    comp=nvl;
    nobj=size(Xcal,1);
    
    % Obter os Escores na validação
    T_val=Xnew_ce/P';
    [Xrecuperado_val]=T_val*P';
    [Xrecuperado_val]=centrar_inverso(Xrecuperado_val,mX);
    Qcont_val=Xnew-Xrecuperado_val; % Residuo de X
      
  
   %%
   
%    Calculo de T2 CALIBRAÇÃO
fvar = sqrt(1./(diag(T'*T)/(size(Xcal,1) - 1)));
Thot = sum((T*diag(fvar)).^2,2);
Tcont = (T*diag(fvar)*P');

%    Calculo de T2 VALIDAÇÃO
Thot_val = sum((T_val*diag(fvar)).^2,2);
Tcont_val = (T_val*diag(fvar)*P');


    % Calculo de T2 limite (Calibração)
lev_conf = 0.95;
if license('test','statistics_toolbox')
    F = finv(lev_conf,comp,nobj-comp);
    tlim = (comp*(nobj - 1)/(nobj - comp))*F;
else
    disp(['Digite o Valor de F com (Namostras-Nº variaveis latentes -1) graus de liberdade ','. ']);
    tlim=input('');   
end
   

%% Calculo de Q Calibração
    for i=1:size(T,1)
    Qres(i) = Qcont(i,:)*Qcont(i,:)';
    end

% Calculo de Q limite (Calibração)

[m,n] = size(Xcal);
s = svd(Qcont');
s = s.^2/(m-1);
Qcont=s;

[m,n] = size(Qcont);

if n>m
  s   = s';
  m   = n;
end

Qcont=s;
t1 = sum(Qcont(0+1:m,1));
t2 = sum(Qcont(0+1:m,1).^2);
t3 = sum(Qcont(0+1:m,1).^3);

ho = 1-2*t1*t3/3/(t2.^2);
if ho<0.001; ho = 0.001; end
ca = sqrt(2)*erfinv(2*lev_conf-1);
term1    = ca*sqrt(2*t2*ho.^2)/t1;
term2    = t2*ho*(ho-1)/(t1.^2);
qlim = t1*(1+term1+term2).^(1/ho);

%% Calculo de Q Validação

    for i=1:size(T_val,1)
    Qres_val(i) = Qcont_val(i,:)*Qcont_val(i,:)';
    end


%% Modelo
model.T = T;
model.T_val=T_val;
model.P = P;
model.U = U;
model.b = B;
model.qlim=qlim;
model.tlim=tlim;
model.Qres=Qres;
model.Qres_val=Qres_val;
model.Tcont_val=Tcont_val;
model.Tcont=Tcont;
model.Thot_val=Thot_val;
model.Thot=Thot
model.Ypre_cal=Yhat;
model.Yprev_val=Ynew_hat;

% Plotando para cada classe na calibração
figure;
title ('Conjunto de Calibração')
for i = 1:size(Ycal,2)
classes{i}=find(Ycal(:,i)==1)
plot(model.Thot(classes{i}),model.Qres(classes{i}),'o')
hold on
end

% figure;plot(model.Thot,model.Qres,'o')
% hold on
vline(model.tlim)
hline (model.qlim)
legend ('show')

% Plotando para cada classe na validação


figure;
title ('Conjunto de Validação')
plot(model.Thot_val,model.Qres_val,'o')
vline(model.tlim)
hline (model.qlim)
legend ('Amostras da Validação')


%% Caso queira plotar a calibração e validação juntas
% NESTE CASO A ULTIMA CLASSE DE AMOSTRAS SERÀ O CONJUNTO DE VALIDAÇÃO
% INTEIRO
% Plotando para cada classe na calibração
figure;
title ('Conjunto de Calibração e Validação')
for i = 1:size(Ycal,2)
classes{i}=find(Ycal(:,i)==1)
plot(model.Thot(classes{i}),model.Qres(classes{i}),'o')
hold on
end
% figure;plot(model.Thot,model.Qres,'o')
% hold on
vline(model.tlim)
hline (model.qlim)
hold on
plot(model.Thot_val,model.Qres_val,'o')
legend ('show')
end
