function [Xp,Xtp]=pretrat(X,Xt,method)
    % rotina para preprocessar dados;
    % rotinas auxiliares: savgol; spdiags.
    % input:
    % X    : matriz de referência para preprocessamento;
    % Xt   : matriz de teste;
    % method : métodos para  preprodessar. É uma célula com os vários métodos de preprocessamento
    %        Ex: method={'center';'auto';'snv';'msc','pareto';'minmax';'minmax2';'norm';'deriv',[7,2,1]}
    %
    % output:
    % Xp   : matriz X preprocessada;
    % Xtp  : matriz Xt preprocessada;
    %
    % Exemplo: method={'msc';'deriv';[7,2,1];'minmax';}
    %          [Xp,Xtp]=pretrat(Xcal,Xprev,method);
    %
    % Paulo R. Filgueiras   -  27/02/2014
    %

[n1,n2]=size(method);
if n2>1; method=method'; end
npret=size(method,1);  % Número de preprocessamento    
for ki=1:npret
    if strcmp(method{ki,1},'center');[X,Xt]=center(X,Xt);
    elseif strcmp(method{ki,1},'auto');[X,Xt]=auto(X,Xt);    
    elseif strcmp(method{ki,1},'snv');[X,Xt]=snv(X,Xt);
    elseif strcmp(method{ki,1},'msc');[X,Xt]=msc(X,Xt);
    elseif strcmp(method{ki,1},'pareto');[X,Xt]=pareto(X,Xt);
    elseif strcmp(method{ki,1},'minmax');[X,Xt]=minmax(X,Xt);
    elseif strcmp(method{ki,1},'minmax2');[X,Xt]=minmax2(X,Xt);
    elseif strcmp(method{ki,1},'deriv');[X,Xt]=derivada(X,Xt,method{ki+1,1});   
    elseif strcmp(method{ki,1},'norm');[X,Xt]=normalizar(X,Xt); 
    elseif strcmp(method{ki,1},'emsc');[X,Xt]=emsc(X,Xt,method{ki+1,1}); 
    elseif strcmp(method{ki,1},'fir');[X,Xt]=fir(X,Xt,method{ki+1,1}); 
    elseif strcmp(method{ki,1},'osc');[X,Xt]=osc(X,Xt,method{ki+1,1},method{ki+2,1}); 
    else ; X=X;Xt=Xt;    
    end
end    
Xp=X;Xtp=Xt;  
    

function [X1,X1p] = center(Xe,Xet)
    % centrar na média
    % input: Xe : matriz X para centrar na média;
    %        Xet: matriz X de teste;
    % output: X1 : matriz Xe centrada na média;
    %         X1p: matriz Xet centrada na média;
    X1=Xe-ones(size(Xe,1),1)*mean(Xe);
    if size(Xe,2)==size(Xet,2)
        X1p=Xet-ones(size(Xet,1),1)*mean(Xe);
    else
        X1p=Xet;
    end
    
function [X1,X1p] = auto(Xe,Xet)
    % autoescalar
    % input: Xe : matriz X para autoescalar;
    %        Xet: matriz X de teste;
    % output: X1 : matriz Xe autoescalada;
    %         X1p: matriz Xet autoescalada;
    X1=(Xe-ones(size(Xe,1),1)*mean(Xe))./(ones(size(Xe,1),1)*std(Xe));
    if size(Xe,2)==size(Xet,2)
        X1p=(Xet-ones(size(Xet,1),1)*mean(Xe))./(ones(size(Xet,1),1)*std(Xe));
    else
        X1p=Xet;
    end
    
function [X1,X1p] = snv(Xe,Xet)
    % suavização normal padrão
    % input: Xe : matriz X para suavizar;
    %        Xet: matriz X de teste;
    % output: X1 : matriz Xe suavizada;
    %         X1p: matriz Xet suavizada;
    X1=(Xe-mean(Xe')'*ones(1,size(Xe,2)))./(std(Xe')'*ones(1,size(Xe,2)));
    if size(Xe,2)==size(Xet,2)
        X1p=(Xet-mean(Xet')'*ones(1,size(Xet,2)))./(std(Xet')'*ones(1,size(Xet,2)));
    else
        X1p=Xet;
    end
    
function [X1,X1p] = msc(Xe,Xet)
    % correção de espalhamento multiplicativo
    % input: Xe : matriz X para correçao;
    %        Xet: matriz X de teste;
    % output: X1 : matriz Xe corrigida;
    %         X1p: matriz Xet corrigida;
    para1=mean(Xe);
    for k2=1:size(Xe,1)
        coef=polyfit(para1,Xe(k2,:),1);
        X1(k2,:)=(Xe(k2,:)-coef(2))/coef(1);
    end
    % corrigindo a matriz de teste
    if size(Xe,2)==size(Xet,2)
    for k3=1:size(Xet,1)
        coef2=polyfit(para1,Xet(k3,:),1);
        X1p(k3,:)=(Xet(k3,:)-coef2(2))/coef2(1);
    end
    else
        X1p=Xet;
    end

function [X1,X1p] = pareto(Xe,Xet)
    % Preprocessamento pareto
    % input: Xe : matriz X para ser preprocessada;
    %        Xet: matriz X de teste;
    % output: X1 : matriz Xe preprocessada;
    %         X1p: matriz Xet preprocessada;
    X1=(Xe-ones(size(Xe,1),1)*mean(Xe))./(ones(size(Xe,1),1)*sqrt(std(Xe)));
    if size(Xe,2)==size(Xet,2)
        X1p=(Xet-ones(size(Xet,1),1)*mean(Xe))./(ones(size(Xet,1),1)*sqrt(std(Xe)));
    else
        X1p=Xet;
    end

function [X1,X1p] = minmax(Xe,Xet)
    % Preprocessamento no intervalo [0,1] de cada variável
    % input: Xe : matriz X para ser preprocessada;
    %        Xet: matriz X de teste;
    % output: X1 : matriz Xe preprocessada;
    %         X1p: matriz Xet preprocessada; 
    X1=(Xe-ones(size(Xe,1),1)*min(Xe))./(ones(size(Xe,1),1)*(max(Xe)-min(Xe)));
    if size(Xe,2)==size(Xet,2)
        X1p=(Xet-ones(size(Xet,1),1)*min(Xe))./(ones(size(Xet,1),1)*(max(Xe)-min(Xe)));
    else
        X1p=Xet;
    end
    
function [X1,X1p] = minmax2(Xe,Xet)
    % Preprocessamento no intervalo [-1, +1] de cada variável
    % input: Xe : matriz X para ser preprocessada;
    %        Xet: matriz X de teste;
    % output: X1 : matriz Xe preprocessada;
    %         X1p: matriz Xet preprocessada; 
    minv=min(Xe);maxv=max(Xe);para1=0.5*(maxv+minv);para2=0.5*(maxv-minv);
    for ki=1:size(Xe,2); X1(:,ki)=(Xe(:,ki)-para1(ki))/para2(ki); end
    if size(Xe,2)==size(Xet,2)
        for ki=1:size(Xet,2); X1p(:,ki)=(Xet(:,ki)-para1(ki))/para2(ki); end
    else
        X1p=Xet;
    end
    
function [X1,X1p] = derivada(Xe,Xet,para1)
    % Derivada com filtro de suavização Savitsky-Golay. 
    % input: Xe : matriz X para ser preprocessada;
    %        Xet: matriz X de teste;
    %      para1: vetor com [janela, grau do polinômio, ordem da derivada]
    % output: X1 : matriz Xe suavizada e diferenciada;
    %         X1p: matriz Xet suavizada e diferenciada;    
    X1=savgol(Xe,para1(1),para1(2),para1(3));   
    if size(Xe,2)==size(Xet,2)
        X1p=savgol(Xet,para1(1),para1(2),para1(3));
    else
        X1p=Xet;
    end
    
function [X1,X1p] = normalizar(Xe,Xet)
    % normalizar
    % input: Xe : matriz X para normalizar;
    %        Xet: matriz X de teste;
    % output: X1 : matriz Xe normalizada;
    %         X1p: matriz Xet normalizada;
    for k2=1:size(Xe,1);   X1(k2,:)=Xe(k2,:)/norm(Xe(k2,:));    end
    
    if size(Xe,2)==size(Xet,2)
        for k2=1:size(Xet,1);   X1p(k2,:)=Xet(k2,:)/norm(Xet(k2,:));    end
    else
        X1p=Xet;
    end
    
    
function [X1,X1p] = emsc(Xe,Xet,opt)
    % Extended Multiplicative Scatter Correction (EMSC) 
    % input: Xe : matriz X para ser preprocessada;
    %        Xet: matriz X de teste;
    %        opt: options
    % output: X1 : matriz Xe suavizada e diferenciada;
    %         X1p: matriz Xet suavizada e diferenciada;    
    para1=mean(Xe);
    if size(Xe,2)==size(Xet,2)
        X1=emscorr(Xe,mean(Xe),opt);
        X1p=emscorr(Xet,mean(Xe),opt);
    else
        X1=emscorr(Xe,mean(Xe),opt);
        X1p=Xet;
    end
    
    
function [X1,X1p] = fir(Xe,Xet,opt)
    % Standardization based on FIR modelling.
    % input: Xe : matriz X para ser preprocessada;
    %        Xet: matriz X de teste;
    %        opt: options
    % output: X1 : matriz Xe suavizada e diferenciada;
    %         X1p: matriz Xet suavizada e diferenciada;    
    if size(Xe,2)==size(Xet,2)
        X1=stdfir(Xe,mean(Xe),opt,1);
        X1p=stdfir(Xet,mean(Xe),opt,1);
    else
        X1=stdfir(Xe,mean(Xe),opt,1);
        X1p=Xet;
    end
    
    
function [X1,X1p] = osc(Xe,Xet,opt1,opt2)
    % OSC
    % input: Xe : matriz X para ser preprocessada;
    %        Xet: matriz X de teste;
    %        opt: options
    % output: X1 : matriz Xe suavizada e diferenciada;
    %         X1p: matriz Xet suavizada e diferenciada;    
    if size(Xe,2)==size(Xet,2)
        [nx,nw,np,nt] = osccalc(Xe,opt1,opt2);
        X1 = oscapp(Xe,nw,np);
        X1p = oscapp(Xet,nw,np);
    else
        [nx,nw,np,nt] = osccalc(Xe,opt1,opt2);
        X1 = oscapp(Xe,nw,np);
        X1p=Xet;
    end
    
    
function [y_hat,D]= savgol(y,width,order,deriv)
%SAVGOL Savitsky-Golay smoothing and differentiation.
%  Inputs are the matrix of ROW vectors to be smoothed (y),
%  and the optional variables specifying the number of points in
%  filter (width), the order of the polynomial (order), and the
%  derivative (deriv). The output is the matrix of smoothed
%  and differentiated ROW vectors (y_hat) and the matrix of 
%  coefficients (cm) which can be used to create a new smoothed/
%  differentiated matrix, i.e. y_hat = y*cm. If number of points,
%  polynomial order and derivative are not specified,
%  they are set to 15, 2 and 0, respectively.
%
%  Example: if y is a 5 by 100 matrix then savgol(y,11,3,1)
%  gives the 5 by 100 matrix of first-derivative row vectors
%  resulting from a 11-point cubic Savitzky-Golay smooth of
%  each row of y.
%
%I/O format is: [y_hat,cm] = savgol(y,width,order,deriv);
%
%See also: MSCORR, SAVGOLCV, SGDEMO, STDFIR, BASELINE, LAMSEL


% Sijmen de Jong Unilever Research Laboratorium Vlaardingen Feb 1993
% Modified by Barry M. Wise 5/94
%         ***   Further modified, 1998-03, Martin Andersson
%         ***   Adjusting the calcn. of the bulk data.
%         ***   Based on calcn. of a sparse derivative matrix (D)

[m,n] = size(y);
y_hat = y;
% set default values: 15-point quadratic smooth
if nargin<4
  deriv= 0;
  disp('  '), disp('Derivative set to zero') 
end
if nargin<3
  order= 2; 
  disp('  '), disp('Polynomial order set to 2')
end
if nargin<2
  width=min(15,floor(n/2)); 
  s = sprintf('Width set to %g',width);
  disp('  '), disp(s)  
end
% In case of input error(s) set to reasonable values
w = max( 3, 1+2*round((width-1)/2) );
if w ~= width
  s = sprintf('Width changed to %g',w);
  disp('  '), disp('Width must be >= 3 and odd'), disp(s)
end
o = min([max(0,round(order)),5,w-1]);
if o ~= order
  s = sprintf('Order changed to %g',o); disp('  ')
  disp('Order must be <= width -1 and <= 5'), disp(s)
end
d = min(max(0,round(deriv)),o);
if d ~= deriv
  s = sprintf('Derivative changed to %g',d); disp('  ')
  disp('Deriviative must be <= order'), disp(s)
end
p = (w-1)/2;
% Calculate design matrix and pseudo inverse
x = ((-p:p)'*ones(1,1+o)).^(ones(size(1:w))'*(0:o));
weights = x\eye(w);
% Smoothing and derivative for bulk of the data
coeff=prod(ones(d,1)*[1:o+1-d]+[0:d-1]'*ones(1,o+1-d,1),1);
D=spdiags(ones(n,1)*weights(d+1,:)*coeff(1),p:-1:-p,n,n);
% Smoothing and derivative for tails 
w1=diag(coeff)*weights(d+1:o+1,:);
D(1:w,1:p+1)=[x(1:p+1,1:1+o-d)*w1]'; 
D(n-w+1:n,n-p:n)=[x(p+1:w,1:1+o-d)*w1]';
% Operate on y using the filtering/derivative matrix, D
y_hat=y*D;

    
    