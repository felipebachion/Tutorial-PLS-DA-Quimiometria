function [Ynew_hat]=previsto_pls(Xcal,Ycal,Xnew,standard,nvl);
% INPUT =
% Xcal = matriz calibração (xcal)
% Ycal = Vetor de calibração (ycal)
% Xnew =  matriz que será prevista Ex: xcal, xval
% standard = 1 para autoescalar os dados e 0 para centrar os dados na média
% nvl = número de variaveis latentes

% OUTPUT = Ynew_hat > valores previstos de Y para o conjunto Xnew

% função para construir o modelo PLS para um dado número de variavies
% latentes e obter os resultados de previsão de Y para o conjunto Xnew
% Caso queira os resultados de calibração utilizo o xcal novamente, exemplo
% para dados centrados na media e número de variaveis latentes = 5
% [ycal_pred]=plspred(xcal,ycal,xcal,0,5);

[n,px]=size(Xcal);   					           
[n,m]=size(Ycal); 

if standard==1						% Caso seja selecionando o autoscalamento
	[Xauto,mX,sX,Xnew_auto]=scale(Xcal,1,Xnew);
	[Yce,mY]=centrar(Ycal,1);				
	[B]=simpls(Xauto,Yce,nvl,Xauto'*Yce,[]);	
	Ynew_hat=Xnew_auto*B;				
	[Ynew_hat]=centrar_inverso(Ynew_hat,mY);
end

if standard==0						
	[Xce,mX,Xnew_ce]=centrar(Xcal,1,Xnew);		
	[Yce,mY]=centrar(Ycal,1);
	[B]=simpls(Xce,Yce,nvl,Xce'*Yce,[]);	        
	Ynew_hat=Xnew_ce*B;			
	[Ynew_hat]=centrar_inverso(Ynew_hat,mY);
end
end
