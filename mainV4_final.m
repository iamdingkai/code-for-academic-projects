

clear;clc;close all;

set(0,'DefaultFigureVisible','on')
%{
cd 'C:\Users\xu8294\Dropbox\Filippo_Kai_Agri_Project(2018)\Computation\Matlab\2018.08.29 recalibrated parameters v3'
addpath 'C:\Users\xu8294\Dropbox\Filippo_Kai_Agri_Project(2018)\Data'; % containing data
addpath 'C:\Users\xu8294\Dropbox\Filippo_Kai_Agri_Project(2018)\Data\CPS_Ag_NonAg_Wage_Gap';
%}

% % %
%   cd '/Users/xu8294/Dropbox/Filippo_Kai_Agri_Project(2018)/Computation/Matlab/2018.08.29 recalibrated parameters v3'
%   addpath '/Users/xu8294/Dropbox/Filippo_Kai_Agri_Project(2018)/Data'; % containing data
%   addpath '/Users/xu8294/Dropbox/Filippo_Kai_Agri_Project(2018)/Data/CPS_Ag_NonAg_Wage_Gap';
% % %}

 cd '/Users/uy3294/Dropbox/Filippo_Kai_Agri_Project(2018)/Computation/Matlab/2018.08.29 recalibrated parameters v3'
 addpath '/Users/uy3294/Dropbox/Filippo_Kai_Agri_Project(2018)/Data'; % containing data
 addpath '/Users/uy3294/Dropbox/Filippo_Kai_Agri_Project(2018)/Data/CPS_Ag_NonAg_Wage_Gap';


set(0,'defaultlinelinewidth',1.3)
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultTextFontSize',12)
set(0,'defaultAxesFontName', 'Times New Roman')
set(0,'defaultTextFontName', 'Times New Roman')
set(0,'DefaultAxesTitleFontWeight','normal');



global pBase_Guess;
pBase_Guess=0.5;


%% Empirical data

% CaliStartY=1986; % 1986 is the first year PSE data is available
% CaliStartY=1996; % 1996 is the first year average ag income data is available


CaliStartY=2001; % starting year for calibration 
%(Ag Emp Share dropped sharply in 1999 and 2000, 
% there is nothing in our model that delivers that. 
% Thefore, we pick 2001 as calibration starting year.)
CaliEndY=2017; % end year for calibration (last year actual data exists for PSE and TFP)
ProjEndY=2020; % end year for projection (from CaliEndY+1 to ProjEndY, we project PSE and TFP series)
T=ProjEndY-CaliStartY+1; % number of periods

% base year for calibration: A=1 for base year, and rho is calibrated to match base year ag-nonag wage ratio 
CaliBaseY=2006; % 2006 is right before financial crisis, so that employment share is not affected by business cycle

Y=(CaliStartY:1:ProjEndY)'; % years including actual data and projection
YA=(CaliStartY:1:CaliEndY)'; % years where actual data points exist
CaliBaseYidx=find(Y==CaliBaseY,1,'first');


% Agri Emp Share
DataAgEmp(:,1)=xlsread('Ag Emp Share US World Bank.xlsx','USA','D1:AD1')';
DataAgEmp(:,2)=xlsread('Ag Emp Share US World Bank.xlsx','USA','D2:AD2')';
Utemp2=DataAgEmp(:,2);
DataAgEmpT=Utemp2((DataAgEmp(:,1)>=CaliStartY)&(DataAgEmp(:,1)<=CaliEndY));

% TFP
DataA=xlsread('Multiple Factor Productivity MFP TFP.xlsx','BLS Data Series','A13:B43');
Utemp2=DataA(:,2);
AT=Utemp2((DataA(:,1)>=CaliStartY)&(DataA(:,1)<=CaliEndY));
AT=AT/AT(CaliBaseYidx);
ATpolfit=polyfit((CaliStartY:1:CaliEndY)',log(AT),1);
ATP=exp(polyval(ATpolfit,Y));
ATP(Y<=CaliEndY)=AT;
AT=ATP;
figure(1);clf;plot(Y(Y<=CaliEndY),AT(Y<=CaliEndY),'bo-',Y,AT,'rx--');

% PSE
DataPSE(:,1)=xlsread('PSE_breakdown_by_commodity.xlsx','Paper','C1:AH1')';
DataPSE(:,2)=xlsread('PSE_breakdown_by_commodity.xlsx','Paper','C8:AH8')';
DataPSE(:,3)=xlsread('PSE_breakdown_by_commodity.xlsx','Paper','C11:AH11')';
Utemp2=DataPSE(:,2);
sT=Utemp2((DataPSE(:,1)>=CaliStartY)&(DataPSE(:,1)<=CaliEndY));
sTpolfit=polyfit((CaliStartY:1:CaliEndY)',log(sT),1);
sTP=exp(polyval(sTpolfit,Y));
sTP(Y<=CaliEndY)=sT;
sT=sTP;
figure(2);clf;plot(Y(Y<=CaliEndY),sT(Y<=CaliEndY),'bo-',Y,sT,'rx--');
Utemp3=DataPSE(:,3);
PSET=Utemp3((DataPSE(:,1)>=CaliStartY)&(DataPSE(:,1)<=CaliEndY));

% Agri Relative Income
DataAgRelInc(:,1)=xlsread('Ag_NonAg_Wage_Gap.xlsx','Kai','A2:A23');
DataAgRelInc(:,2)=xlsread('Ag_NonAg_Wage_Gap.xlsx','Kai','B2:B23');
Utemp2=DataAgRelInc(:,2);
DataAgRelIncT=Utemp2((DataAgRelInc(:,1)>=CaliStartY)&(DataAgRelInc(:,1)<=CaliEndY));
AvgWageRatioBase=Utemp2(DataAgRelInc(:,1)==CaliBaseY); % Lagakos Waugh uses this target

% Agri Relative Price (using relative agriculture deflator)
DataP(:,1)=xlsread('ag_deflator_ag_import_ag_export_ag_emp.xlsx','Sheet1','A2:A32');
DataP(:,2)=xlsread('ag_deflator_ag_import_ag_export_ag_emp.xlsx','Sheet1','D2:D32');
temp=DataP(:,2);
DataPT=temp((DataP(:,1)>=CaliStartY)&(DataP(:,1)<=CaliEndY));
DataPT=[DataPT;nan;nan;nan;nan;]; % relative ag deflator only goes to 2016 while Y goes to 2020


%% parameters taken from Lagakos-Waugh (thetaa and thetan)
% parameters for productivity draws
thetaa=5.3;
thetan=2.7;


%% calibrating "rho" to match "current ag emp share and current ag non-ag income ratio"
rhoLagakos=3.5; % from Lagakos and Waugh, as initial guess
AgEmpBase=DataAgEmpT(CaliBaseYidx)/100; % ag emp share in base year
sBase=sT(CaliBaseYidx); % PSE in base year
ABase=1; % normalization
N=1000; % how many grids to use to discretize za and zn

% used for plotting utility and utility changes
Frac=1; % calculating/plotting utilities of 1:1:N/Frac
D1=1:1:N/Frac;
D2=1:1:N/Frac;

% D1=1:Frac:N;
% D2=1:Frac:N;

rho=fzero(@(rho)Fun_Cali_rho(rho,AgEmpBase,sBase,AvgWageRatioBase,thetaa,thetan,N),...
    [rhoLagakos*0.7,rhoLagakos*1.5],optimset('Display','iter'));
[~,pBase,za1,zn1,za2,zn2,g2,zaV,znV,gV]...
    =Fun_Cali_rho(rho,AgEmpBase,sBase,AvgWageRatioBase,thetaa,thetan,N);

%% estimated net export as a function of domestic relative price of ag
gam1=2.702; % elasticity of "NX/ Ag GDP" to "Ag Deflator/Deflator", gam1 is estimated from data
gam0=fsolve(@(gam0)gam0*pBase^(-gam1)-0.15,0.15,optimset('Display','iter')); % gam0 is calibrated to match Base year NX to Ag GDP ratio of 15%

%gam1=0;
%gam0=0;


%% calibrating "nu" to match "LR agri emp share"
AgEmpLR=0.5/100;    % Lagakos Waugh uses this target of long run ag emp share
sLR=0;                          % long run level of subsidy, assume is 0
AgEmpLR_handle=@(p)sum(gV.*(p*(1+sLR).*zaV>=znV)); % percentage emp in agri
pLR=fzero(@(p)100*(AgEmpLR_handle(p)-AgEmpLR),[0 1],optimset('Display','iter')); % calculate the relative ag price that matches long-run ag emp share
yaLR=sum(gV.*(pLR.*(1+sLR).*zaV>=znV).*zaV); % long run agri production normalized by A
ynLR=sum(gV.*(pLR.*(1+sLR).*zaV<znV).*znV); % long run non-agri production normalized by A
caLR=yaLR*(1-gam0*pLR^(-gam1)); % long run agri consumption normalized by A
nu=ynLR/caLR/pLR; % calculate the nu that is consistent with HH first order condition

%% calibrating "a_bar" to match "base year consumer FOC"
yaBase=sum(gV.*(pBase.*(1+sBase).*zaV>=znV).*zaV); % base year agri production normalized by A
ynBase=sum(gV.*znV.*(pBase*(1+sBase).*zaV<znV)); % base year non-agri production normalized by A
caBase=yaBase*(1-gam0*pBase^(-gam1)); % base year agri consumption normalized by A
abar=ABase*(caBase-ynBase/(pBase*nu)); % calculate the abar consistent with Base year FOC

% utility function 
u=@(ca,cn)log(ca-abar)+nu*log(cn);
u2=@(p,y)log((y-p*abar)./p./(1+nu))+nu*log(nu*(y-p*abar)./(1+nu)); % u2(p,y), y is income after subsidy and tax



%% Benchmark
% Calculate the time series of consumption, employment, GDP from starting year to ending year
[B.yaT,B.ynT,B.AgEmpT,B.NAgEmpT,B.tauT]=deal(nan(T,1));
B.AgriIndV2=nan(N*N,T);
B.IncT3=nan(N,N,T); % bring home income
B.pT=5*ones(T,1); % initial guess
pLB=1e-3;
pUB=pBase*1.5;

for t=1:T
    s=sT(t);
    A=AT(t);
    
    yn=@(p)A*sum((p.*(1+s).*zaV<znV).*znV.*gV); % non agri output
    ya=@(p)A*sum((p.*(1+s).*zaV>znV).*zaV.*gV); % agri output
    ca=@(p)ya(p).*(1-gam0*p.^(-gam1)); % domestic consumption of agri product
    
    p=fzero(@(p)yn(p)-nu*p*(ca(p)-abar),[pLB pUB],optimset('display','iter','tolX',pLR/1e3)); % p solves FOC of households
    pLB=max(0,p*0.5); % new lower bound
    pUB=p*1.3; % new upper bound
    
    B.AgriIndV=(p*(1+s)*zaV>=znV);                   % indicator for working in Ag
    B.yaT(t)=A*sum(B.AgriIndV.*zaV.*gV);             % agri output in units of crops
    B.ynT(t)=A*sum((1-B.AgriIndV).*znV.*gV);       % non-agri output in units of cars
    B.pT(t)=p;                                                                 % relative price of Ag goods in year t
    
    B.AgEmpT(t)=sum(B.AgriIndV.*gV);                  % total/proportion of employment in agri in period t
    B.NAgEmpT(t)=sum((1-B.AgriIndV).*gV);         % total/proportion of employment in non-agri in period t
    
    B.tauT(t)=s*p*B.yaT(t)/((1+s)*p*B.yaT(t)+B.ynT(t));  % proportional tax rate on everyone to finance PSE spending
    B.IncT3(:,:,t)=A*(1-B.tauT(t))*max(zn2,(1+s)*p*za2); % income after tax and subsidy (in units of cars)
end

B.UT3=u2(repmat(permute(B.pT,[2 3 1]),[N N 1]),B.IncT3);                                      % utility of individual (za,zn) in period t
B.AvgAgInc_Include_Subsidy_T=B.pT.*B.yaT.*(1+sT).*(1-B.tauT)./B.AgEmpT;    % average agri income after subsidy and tax in units of cars
B.AvgAgInc_Exclude_Subsidy_T=B.pT.*B.yaT./B.AgEmpT;                                     % average agri income after subsidy and tax in units of cars
B.AvgNAgIncT=B.ynT.*(1-B.tauT)./B.NAgEmpT;                                                      % in units of cars    
B.GDPT=B.ynT+B.pT.*B.yaT;


% plot demand and supply curves
pVec=linspace(B.pT(1)/2,B.pT(1)*2,20);
[yaVec,ynVec,caVec]=deal(nan(size(pVec)));
for i=1:20
    p=pVec(i);
    s=sT(1);
    A=AT(1);
    AgriIndV_temp=(p*(1+s)*zaV>=znV);                   % indicator for working in Ag
    yaVec(i)=A*sum(AgriIndV_temp.*zaV.*gV);             % agri output in units of crops
    IncTemp=B.IncT3(:,:,1);
    caV=(IncTemp(:)+p*nu*abar)/(p*(1+nu));
    caVec(i)=sum(caV.*gV); % total domestic demand for agriculture product in period 1
end

figure(1000);clf;
plot(caVec,pVec,'b-',yaVec,pVec,'r',yaVec.*(1-gam0*pVec.^(-gam1)),pVec,'r--');
box on;grid on;legend('D','S');
ylabel P; xlabel Q;
% xlim([0 0.4]);
xlim([0.02 0.031]);
ylim([0.3 0.55]);


figure(1001);clf;
plot(pVec,caVec,'b-',pVec,yaVec,'r',pVec,yaVec.*(1-gam0*pVec.^(-gam1)),'r--');
box on;grid on;legend('D','S');
xlabel P; ylabel Q;


% #################################################
% plots for the draft

% Exogenous Sequence to Feed into Model: PSE
figure(1);clf;
plot(YA,PSET*100,'bo-',Y,sT*100,'rx--');box on;grid on;
xlabel Year;ylabel ('% of Agriculture Value-Added');
legend('PSE: All Components','PSE: Distortionary Components');
% title 'Agriculture Subsidy (PSE)';
xlim([CaliStartY CaliEndY]);
figure_name = 'Distortionary_PSE';
r=150;set(gcf,'PaperPosition',[0 0 600 500]/r);
% print('-depsc2',figure_name);

% Exogenous Sequence to Feed into Model: TFP
figure(2);clf;
subplot(1,2,1);plot(Y,sT*100,'bo-');box on;grid on;
xlabel Year;ylabel ('% of Agriculture Value-Added');
title 'Agriculture Subsidy (PSE)';
xlim([CaliStartY CaliEndY]);
subplot(1,2,2);plot(Y,AT,'bo-');box on;grid on;
xlabel Year;ylabel ('Normalized by 2006 Value');
xlim([CaliStartY CaliEndY]);
title 'Total Factor Productivity';
figure_name = 'Subsidy_TFP';
r=150;set(gcf,'PaperPosition',[0 0 1200 500]/r);
% print('-depsc2',figure_name);

% Agriculture Employment Share
figure(3);clf;
plot(Y,B.AgEmpT*100,'r--',YA,DataAgEmpT,'b-');box on;grid on;hold on;
plot(CaliBaseY,B.AgEmpT(Y==CaliBaseY)*100,'ko','MarkerFaceColor', 'k');
xlabel Year;ylabel ('Percentage');
xlim([CaliStartY CaliEndY]);
legend('Model','Data');
figure_name = 'Agri_Emp_Share';
r=150;set(gcf,'PaperPosition',[0 0 600 500]/r);
% print('-depsc2',figure_name);

% Agriculture Relative Income
figure(4);clf;
plot(Y,B.AvgAgInc_Include_Subsidy_T./B.AvgNAgIncT*100,'r--');box on;grid on;hold on;
plot(YA,DataAgRelIncT*100,'b-');box on;grid on;
temp=B.AvgAgInc_Include_Subsidy_T(Y==CaliBaseY)./B.AvgNAgIncT(Y==CaliBaseY)*100;
plot(CaliBaseY,temp,'ko','MarkerFaceColor', 'k');
ylim([66 80]);xlim([CaliStartY CaliEndY]);
xlabel Year;ylabel ('Percentage');
legend('Model','Data','location','SE');
figure_name = 'Agri_Relative_Income';
r=150;set(gcf,'PaperPosition',[0 0 600 500]/r);
% print('-depsc2',figure_name);

% Relative Price of Agri Output
figure(5);clf;
plot(Y, B.pT,'ro-');box on;grid on;
hold on;
plot(Y, B.pT.*(1+sT),'b+-');
xlabel Year;ylabel ('Relative Price');
xlim([CaliStartY CaliEndY]);
legend('p','p*(1+s)','location','SE');
figure_name = 'Relative_Price_Agri_Output';
r=150;set(gcf,'PaperPosition',[0 0 600 500]/r);
% print('-depsc2',figure_name);

% utility of everyone (over the entire time frame) as a function of (za,zn)
figure(7);clf;
Period_of_Interest=1;
Utemp2=B.UT3(:,:,Period_of_Interest);
Utemp2(imag(Utemp2)~=0)=nan;
surf(za2(D1,D2),zn2(D1,D2),Utemp2(D1,D2)); % first period
alpha 0.8;
xlabel z^a;ylabel z^n;zlabel Utility;
Cutofftemp2=reshape(B.AgriIndV(:,Period_of_Interest),[N N]); % agriculture occupation indicator
Cutoff_y=nan(N,1);
Cutoff_z=nan(N,1);
for i=1:N % over zn
   idx=find(Cutofftemp2(i,:)==0,1,'first'); % first the last/largest za such that households still stay in non-agri 
   Cutoff_y(i)=zn1(idx);
   Cutoff_z(i)=Utemp2(i,idx);
end
hold on;
h=plot3(za1,Cutoff_y,Cutoff_z,'r-','linewidth',6);
uistack(h,'top');
xlim([0 5]);ylim([0 8]);
figure_name = 'Utility by Productivity Draws';
r=150;set(gcf,'PaperPosition',[0 0 600 500]/r);
% print('-depsc2',figure_name);

% relative ag price: model vs data
figure(8);clf;
Per_Dev=@(x)(x-x(1))./x(1)*100; % percentage deviation from starting value
% plot(Y,DataPT,'bo-',Y,B.pT,'r+-',Y,B.pT.*(1+sT),'mx--');box on;grid on;
plot(Y,Per_Dev(DataPT),'bo-',Y,Per_Dev(B.pT),'r+-',Y,Per_Dev(B.pT.*(1+sT)),'mx--');box on;grid on;
xlabel Year;ylabel ('Relative Price of Agriculture Output');
xlim([CaliStartY CaliEndY]);
legend('Data','Model p','Model p(1+s)','location','NW');
xlim([CaliStartY CaliEndY]);
figure_name = 'Data_p_Model_p';
r=150;set(gcf,'PaperPosition',[0 0 600 500]/r);
% print('-depsc2',figure_name);

% end plots for the draft
% #################################################


%% Counterfactual I: Remove Subsidies, Add Lump-Sum Transfer (LS) to Keep Utility Constant Period by Period
[LS.yaT,LS.ynT,LS.AgEmpT,LS.NAgEmpT,LS.SubAgT]=deal(nan(T,1));
LS.IncT3=nan(N,N,T); % bring home income
LS.pT=5*ones(T,1); % initial guess
pLB=1e-3;
pUB=pBase*2;

for t=1:T
    A=AT(t);
    
    yn=@(p)A*sum((p.*zaV<znV).*znV.*gV); % employment indicator is based on no subsidy
    ya=@(p)A*sum((p.*zaV>znV).*zaV.*gV);
    ca=@(p)ya(p).*(1-gam0*p.^(-gam1)); % domestic consumption of agri product
    
    p=fzero(@(p)yn(p)-nu*p*(ca(p)-abar),[pLB pUB],optimset('display','iter','tolX',pLR/1e3));
    pLB=max(0,p*0.9);
    pUB=p*1.1;
    
    LS.AgriInd2=(p*za2>zn2); % agriculture indicator
    LS.AgriIndV=(p*zaV>znV); % agriculture indicator
    LS.yaT(t)=A*sum(LS.AgriIndV.*zaV.*gV);          % in units of crops
    LS.ynT(t)=A*sum((1-LS.AgriIndV).*znV.*gV);      % in units of cars
    LS.pT(t)=p; 
    
    LS.AgEmpT(t)=sum(LS.AgriIndV.*gV);
    LS.NAgEmpT(t)=sum((1-LS.AgriIndV).*gV);
    
    LS.IncT3(:,:,t)=p*abar+exp((B.UT3(:,:,t)+log(p)+(1+nu)*log(1+nu)-nu*log(nu))/(1+nu));  % bring-home income needed to generate utility level as before
    LS.SubT3(:,:,t)=LS.IncT3(:,:,t)-A*max(zn2,za2.*p); % required lump-sum transfer to individual (za,zn) in period t
    
    Utemp2=LS.SubT3(:,:,t);
    
    LS.SubAgT(t)=sum(Utemp2(:).*LS.AgriIndV.*gV); % total lump-sum transfer to farmers paid out by the planner
    LS.RevT(t)=-sum(sum(LS.SubT3(:,:,t).*g2)); % total revenue to gov't over time
    
    %LS.AvgAgIncT(t)=LS.pT(t)*LS.yaT(t)/LS.AgEmpT(t); % in units of cars
    %LS.AvgNAgIncT(t)=LS.ynT(t)/LS.NAgEmpT(t); % in units of cars    
    
    temp=LS.IncT3(:,:,t);
    LS.AvgAgIncT(t)=sum(temp(:).*LS.AgriIndV.*gV)./LS.AgEmpT(t); % in units of cars
    LS.AvgNAgIncT(t)=sum(temp(:).*(1-LS.AgriIndV).*gV)./(LS.NAgEmpT(t));

    
    LS.GDPT(t)=LS.ynT(t)+LS.pT(t)*LS.yaT(t);
    LS.AgGDPT(t)=LS.pT(t)*LS.yaT(t); % Agri GDP in units of cars
end




% #################################################
% plots for the draft

% Relative Price of Agri Output (before and after subsidy removal)
figure(21);clf;
plot(Y,B.pT,'magenta-.',Y,B.pT.*(1+sT),'b-',Y,LS.pT,'r--');
box on;grid on;
text(2002.5,0.527,'p^a(1+s)','fontsize',14);

text(2002.5,0.434,'{p^a}_{Partial Equilibrium}','fontsize',14);
% text(2002.5,0.428,'(Partial Equilibrium)','fontsize',13);

text(2002.5,0.503,'{p^a}_{General Equilibrium}','fontsize',14);
% text(2002.5,0.497,'General Equilibrium','fontsize',13);

hold on;
% PE=quiver(2000,B.pT(Y==2002).*(1+sT(Y==2002)),0,-B.pT(Y==2002).*sT(Y==2002),'color','r','linewidth',2);
% legend('PSE: p','PSE: p(1+s)','Replace PSE with LS Transfer: p','location','SE');
xlabel Year;ylabel ('Relative Price of Agriculture Output');
xlim([CaliStartY CaliEndY]);
figure_name = 'Relative_Price_Agri_Output_GE';
r=150;set(gcf,'PaperPosition',[0 0 600 500]/r);
% print('-depsc2',figure_name);

% Relative Price of Agri Output (before and after subsidy removal) ####OLD
figure(221);clf;
Norm=@(x)(x-B.pT.*(1+sT))./(B.pT.*(1+sT))*100;
plot(Y,Norm(B.pT),'b+-',Y,Norm(B.pT.*(1+sT)),'bo-',Y,Norm(LS.pT),'rv-');
legend('PSE: p','PSE: p(1+s)','Replace PSE with LS Transfer: p','location','SE');box on;grid on;
xlabel Year;ylabel ('Relative Price of Agriculture Output');
xlim([CaliStartY CaliEndY]);
% figure_name = 'Relative_Price_Agri_Output_GE';
% r=150;set(gcf,'PaperPosition',[0 0 600 500]/r);
% print('-depsc2',figure_name);

% Average Ag Income (before and after subsidy removal)
figure(22);clf; 
plot(Y,B.AvgAgInc_Include_Subsidy_T,'bo-',Y,LS.AvgAgIncT,'r+-');grid on;box on;
legend('PSE','LS Transfer','location','SE');box on;grid on;
xlim([CaliStartY CaliEndY]);
xlabel Year;ylabel ('Average Agriculture Income (Subsidy Included)');
figure_name = 'Average_Agri_Income_GE';
r=150;set(gcf,'PaperPosition',[0 0 600 500]/r);
% print('-depsc2',figure_name);

% Relative Income (before and after subsidy removal)
figure(24);clf;
plot(Y,B.AvgAgInc_Include_Subsidy_T./B.AvgNAgIncT*100,'bo-')
hold on;
plot(Y,LS.AvgAgIncT./LS.AvgNAgIncT*100,'r+-');grid on;box on;
ylim([68 72]);xlim([CaliStartY CaliEndY]);
legend('PSE','LS Transfer','location','SE');box on;grid on;
xlabel Year;ylabel ('Relative Hourly Income');
figure_name = 'Relative_Agri_Income_GE';
r=150;set(gcf,'PaperPosition',[0 0 600 500]/r);
% print('-depsc2',figure_name);

% Agri Emp Share (before and after subsidy removal)
figure(28);clf;
plot(Y,B.AgEmpT*100,'bo-',Y,LS.AgEmpT*100,'r+-');grid on;box on;
legend('PSE','LS Transfer','location','NE');box on;grid on;
xlabel Year;ylabel ('Agriculture Employment Share');
xlim([CaliStartY CaliEndY]);
figure_name = 'Agri_Emp_Share_GE';
r=150;set(gcf,'PaperPosition',[0 0 600 500]/r);
% print('-depsc2',figure_name);

% Subsidy for each type (za,zn)
figure(29);clf;
SubsidyPeriod1=LS.SubT3(:,:,1);
SubsidyPeriod1(imag(SubsidyPeriod1)~=0)=nan;
LS.IncPeriod1=LS.IncT3(:,:,1);
LS.IncPeriod1(imag(LS.IncPeriod1)~=0)=nan;
D1=1:50:N;D2=1:50:N;
surf(za2(D1,D2),zn2(D1,D2),SubsidyPeriod1(D1,D2)./LS.IncPeriod1(D1,D2)*100)
xlabel z^a;ylabel z^n;zlabel 'Subsidy (Percentage of Income)'; 
figure_name = 'Agri_Subsidy_by_za_zn_Distribution';
r=150;set(gcf,'PaperPosition',[0 0 600 500]/r);
% print('-depsc2',figure_name);


% End Plots for the Draft
% #################################################


%% Counterfactual II: Remove Subsidies, Add Lump-Sum Transfer to Keep Utility Constant Period by Period BUT Keep p fixed at the level in the Benchmark (Partial Equilibrium analysis, PE)
[PE.yaT,PE.ynT,PE.AgEmpT,PE.NAgEmpT,PE.SubAgT]=deal(nan(T,1)); 
PE.IncT3=nan(N,N,T); % bring home income 
PE.pT=B.pT; % holding price fixed at the benchmark level

for t=1:T
    A=AT(t);
    p=PE.pT(t); % holding price fixed at the benchmark level
    
    PE.AgriIndV=(p*zaV>znV);                  % indicator of sector choice: holding p fixed and removing s (s=0) should lead to relocation of ag workers to non-ag sector
    PE.yaT(t)=A*sum(PE.AgriIndV.*zaV.*gV);          % Ag output: in units of crops
    PE.ynT(t)=A*sum((1-PE.AgriIndV).*znV.*gV);      % Non-Ag output: in units of cars
    
    PE.AgEmpT(t)=sum(PE.AgriIndV.*gV);              % Ag Emp Share
    PE.NAgEmpT(t)=sum((1-PE.AgriIndV).*gV);         % Non-Ag Emp Share
    
    PE.IncT3(:,:,t)=p*abar+exp((B.UT3(:,:,t)+log(p)+(1+nu)*log(1+nu)-nu*log(nu))/(1+nu));  % bring-home income needed to generate utility level as before
    PE.SubT3(:,:,t)=PE.IncT3(:,:,t)-A*max(zn2,za2.*p); % required lump-sum transfer to individual (za,zn) in period t
    
    Utemp2=PE.SubT3(:,:,t);
    
    PE.SubAgT(t)=sum(Utemp2(:).*PE.AgriIndV.*gV); % total lump-sum transfer paid out by the planner
end
PE.RevT=-squeeze(sum(sum(bsxfun(@times,PE.SubT3,g2),1),2)); % total revenue to gov't over time

PE.AvgAgIncT=PE.pT.*PE.yaT./PE.AgEmpT;  % Avg Ag Income: in units of cars
PE.AvgNAgIncT=PE.ynT./PE.NAgEmpT;       % Avg Non-Ag Income: in units of cars
PE.GDPT=PE.ynT+PE.pT.*PE.yaT;           % GDP: in units of cars
PE.AgGDPT=PE.pT.*PE.yaT;                % Agri GDP: in units of cars


% #################################################
% plots for the draft

% Average Income (before and after subsidy removal)
% Blue: average income of farmers under PSE (including subsidy)
% Purple: remove PSE, say 15%, price remains fixed at old level (Partial
% Equilibrium), so income of each farmer drops by 15%. Which should have 
% lead to 15% drop in average ag income. However, due to selection, some
% farmers leave ag sector, which drives up avg ag income.
figure(26);clf;
plot(Y,B.AvgAgInc_Include_Subsidy_T,'bo-',Y,LS.AvgAgIncT,'r+-',Y,PE.AvgAgIncT,'mx--');grid on;box on;
legend('PSE','LS Transfer','LS Transfer, Price Fixed','location','SE');box on;grid on;
xlabel Year;ylabel ('Average Agriculture Income');
xlim([CaliStartY CaliEndY]);
figure_name = 'Average_Agri_Income_PE';
r=150;set(gcf,'PaperPosition',[0 0 600 500]/r);
% print('-depsc2',figure_name);

% Relative Income (before and after subsidy removal)
figure(27);clf;
plot(Y,B.AvgAgInc_Include_Subsidy_T./B.AvgNAgIncT*100,'bo-');
hold on;
plot(Y,LS.AvgAgIncT./LS.AvgNAgIncT*100,'r+-');
plot(Y,PE.AvgAgIncT./PE.AvgNAgIncT*100,'mx--');
grid on;box on;xlim([CaliStartY CaliEndY]);
legend('PSE','LS Transfer','LS Transfer, Price Fixed','location','NE');box on;grid on;
xlabel Year;ylabel ('Relative Hourly Income');
figure_name = 'Relative_Agri_Income_PE';
r=150;set(gcf,'PaperPosition',[0 0 600 500]/r);
% print('-depsc2',figure_name);

% Agri Emp Share (before and after subsidy removal)
figure(28);clf;
plot(Y,B.AgEmpT*100,'b-',Y,LS.AgEmpT*100,'r--',Y,PE.AgEmpT*100,'m-.');grid on;box on;
% legend('PSE','LS Transfer','LS Transfer, Price Fixed','location','SE');box on;grid on;
xlabel Year;ylabel ('Agriculture Employment Share');
xlim([CaliStartY CaliEndY]);
figure_name = 'Agri_Emp_Share_PE';
r=150;set(gcf,'PaperPosition',[0 0 600 500]/r)

text(2002.5,1.73,'with s_t','fontsize',14);
text(2002.5,1.3,'w/o s_t, General Equilibrium','fontsize',14);
text(2002.5,0.5,'w/o s_t, Partial Equilibrium','fontsize',14);

% print('-depsc2',figure_name);

% Fiscal Savings (before and after subsidy removal)
% Idea: compare Benchmark PSE and Lump-Sum transfer. Due to GE effect,
% distortive effect of PSE on sectoral allocation is small. Therefore, welfare 
% loss is small. Note WELFARE LOSS = FISCAL SAVING (dual). So fiscal 
% saving is small when ag price adjusts. However, distortion is large when
% price is held fixed hence large fiscal saving as well. 
figure(28);clf;
plot(Y,LS.RevT./LS.AgGDPT*100,'bo-',Y,PE.RevT./PE.AgGDPT*100,'r+-',Y,sT*100,'mx--');grid on;box on;
legend('LS Transfer','LS Transfer, Price Fixed','PSE','location','NE');box on;grid on;
xlabel Year;ylabel ('Savings: Percentage of Agri GDP');
xlim([CaliStartY CaliEndY]);
figure_name = 'Savings_Perc_Ag_GDP_PE';
r=150;set(gcf,'PaperPosition',[0 0 600 500]/r);
% print('-depsc2',figure_name);

figure(280);clf;
plot(Y,LS.RevT./LS.AgGDPT*100,'b-',Y,PE.RevT./PE.AgGDPT*100,'r--');grid on;box on;
legend('General Equilibrium','Partial Equilibrium','location','NE');box on;grid on;
xlabel Year;ylabel ('Savings: % of Agricultural GDP');
xlim([CaliStartY CaliEndY]);
figure_name = 'Savings_Perc_Ag_GDP_GE_PE';
r=150;set(gcf,'PaperPosition',[0 0 600 500]/r);
% print('-depsc2',figure_name);


% Fiscal Savings (before and after subsidy removal)
% Idea: compare Benchmark PSE and Lump-Sum transfer. Due to GE effect,
% distortive effect of PSE on sectoral allocation is small. Therefore, welfare 
% loss is small. Note WELFARE LOSS = FISCAL SAVING (dual). So fiscal 
% saving is small when ag price adjusts. However, distortion is large when
% price is held fixed hence large fiscal saving as well. 
figure(29);clf;
plot(Y,LS.RevT./(LS.AgGDPT.*sT')*100,'b-',Y,PE.RevT./(PE.AgGDPT.*sT)*100,'r--');grid on;box on;
legend('General Equilibrium','Partial Equilibrium','location','NE');box on;grid on;
xlabel Year;ylabel ('Savings: % of Current Transfers s_t');
xlim([CaliStartY CaliEndY]);
figure_name = 'Savings_Perc_Ag_PSE_GE_PE_Perc_of_PSE';
r=150;set(gcf,'PaperPosition',[0 0 600 500]/r);
% print('-depsc2',figure_name);



% End Plots for the Draft
% #################################################


%% Counterfactual III: Remove PSE and add no subsidy (ADD NOTHING), look at welfare of everyone (za,zn)
% This is the first counterfactual in the paper - Filippo 11/15/2018
% Argue that PSE is NOT a subsidy to farmers, but a subsidy to the poor
[NO.yaT,NO.ynT,NO.AgEmpT,NO.NAgEmpT,NO.SubAgT]=deal(nan(T,1));
NO.AgriIndV2=nan(N*N,T);    
NO.IncT3=nan(N,N,T); % bring home income
NO.pT=5*ones(T,1); % initial guess
pLB=1e-3;
pUB=pBase*1.5;

for t=1:T
    A=AT(t);
    
    yn=@(p)A*sum((p.*zaV<znV).*znV.*gV);
    ya=@(p)A*sum((p.*zaV>znV).*zaV.*gV);
    ca=@(p)ya(p).*(1-gam0*p.^(-gam1)); % domestic consumption of agri product
    
    p=fzero(@(p)yn(p)-nu*p*(ca(p)-abar),[pLB pUB],optimset('display','iter','tolX',pLR/1e3));
    pLB=max(0,p*0.9);
    pUB=p*1.1;
    
    NO.AgriIndV=(p*zaV>znV);     
    NO.AgriIndV2(:,t)=(p*zaV>znV); 

    NO.yaT(t)=A*sum(NO.AgriIndV.*zaV.*gV); % in units of crops
    NO.ynT(t)=A*sum((1-NO.AgriIndV).*znV.*gV);      % in units of cars
    NO.pT(t)=p;
    
    NO.AgEmpT(t)=sum(NO.AgriIndV.*gV);
    NO.NAgEmpT(t)=sum((1-NO.AgriIndV).*gV);
    NO.IncT3(:,:,t)=A*max(zn2,za2.*p); % income
end
NO.UT3=u2(repmat(permute(NO.pT,[2 3 1]),[N N 1]),NO.IncT3); % utility of individual (za,zn) in period t
NO.AvgAgIncT=NO.pT.*NO.yaT./NO.AgEmpT; % in units of cars
NO.AvgNAgIncT=NO.ynT./NO.NAgEmpT; % in units of cars
NO.GDPT=NO.ynT+NO.pT.*NO.yaT;
NO.AgGDPT=NO.pT.*NO.yaT; % Agri GDP

% Relative Price of Agri Output (before and after subsidy removal)
figure(45);clf;
plot(Y,B.pT,'b+-',Y,B.pT.*(1+sT),'bo-',Y,NO.pT,'rv-');
box on;grid on;
xlim([CaliStartY CaliEndY]);
text(2001,0.53,'PSE: p(1+s_t)');
text(2001,0.43,'PSE: p');
text(2001,0.49,'PSE Removed: p');
hold on;
xlabel Year;ylabel ('Relative Price of Agriculture Output');
figure_name = 'Relative_Price_Agri_Output_PSE_vs_NO';
r=150;set(gcf,'PaperPosition',[0 0 600 500]/r);
% print('-depsc2',figure_name);


% utility of everyone (over the entire time frame) as a function of (za,zn)
figure(41);clf;Period_of_Interest=1;
Utemp2=B.UT3(:,:,Period_of_Interest);
Utemp2(imag(Utemp2)~=0)=nan;
% subplot(1,2,1);
surf(za2(D1,D2),zn2(D1,D2),Utemp2(D1,D2)); % first period
xlabel z^a;ylabel z^n;zlabel Utility;
hold on;
NO.Utemp2=NO.UT3(:,:,Period_of_Interest);
NO.Utemp2(imag(NO.Utemp2)~=0)=nan;
% subplot(1,2,2);
surf(za2(D1,D2),zn2(D1,D2),Utemp2(D1,D2));
%{
% change in utility of everyone (over the entire time frame) as a function of (za,zn)
figure(42);clf;Period_of_Interest=1;
Utemp2=B.UT3(:,:,Period_of_Interest);
Utemp2(imag(Utemp2)~=0)=nan;
NO.Utemp2=NO.UT3(:,:,Period_of_Interest);
NO.Utemp2(imag(NO.Utemp2)~=0)=nan;
Diff2=NO.Utemp2-Utemp2;
surf(za2(D1,D2),zn2(D1,D2),Diff2(D1,D2)); % first period
hold on;
% surf(za2(D1,D2),zn2(D1,D2),zeros(length(D1),length(D2)));alpha 0.5;
xlabel z^a;ylabel z^n;zlabel 'Change in Utility';
figure_name = 'Change_In_Utility_If_Completely_Remove_Subsidy_GE';
r=150;set(gcf,'PaperPosition',[0 0 600 500]/r)
% print('-depsc2',figure_name);
%}

% change in utility of everyone (over the entire time frame) as a function
% of (za,zn) with cutoff line added
figure(43);clf;Period_of_Interest=1;
Utemp2=B.UT3(:,:,Period_of_Interest);
Utemp2(imag(Utemp2)~=0)=nan;
NO.Utemp2=NO.UT3(:,:,Period_of_Interest);
NO.Utemp2(imag(NO.Utemp2)~=0)=nan;
Diff2=NO.Utemp2-Utemp2;
PercDiff2=Diff2./(-(Utemp2-max(Utemp2(:))))*100;  % Utemp2 changes from negative to positive so need to normalize
% h=surf(za2(D1,D2),zn2(D1,D2),PercDiff2(D1,D2)); % first period
h=surf(za2,zn2,PercDiff2); % first period
% colormap polarmap
hold on;
view([-45 40]);
colorbar;
caxis([-10 0]);
xlim([0 za1(N/4)]);
ylim([0 za1(N/4)]);
zlim([-8 0.5]);
set(h,'edgecolor','none');
xlabel z^a;ylabel z^n;zlabel 'Percentage Change in Utility';title('General Equilibrium Effect');
Period_of_Interest=1;
Cutofftemp2=reshape(NO.AgriIndV2(:,Period_of_Interest),[N N]); % agriculture occupation indicator
Cutoff_y=nan(N,1);
Cutoff_z=nan(N,1);
for i=1:N % over zn
   idx=find(Cutofftemp2(i,:)==0,1,'first'); % first the last/largest za such that households still stay in non-agri 
   Cutoff_y(i)=zn1(idx);
   Cutoff_z(i)=PercDiff2(i,idx);
end
ciccio_y=movmean(Cutoff_y,20);
ciccio_z=movmean(Cutoff_z,50);
hold on;
%h=plot3(za1,Cutoff_y,Cutoff_z,'r-','linewidth',0.1);
h=plot3(za1,ciccio_y,ciccio_z,'r','linewidth',1.5);
clear ciccio_y
uistack(h,'top');
Cutoff_0_x=za1;
Cutoff_0_y=nan(N,1);
Cutoff_0_z=nan(N,1);
for i=1:N
    temp=Diff2(i,:);
    idx=find(temp>=0,1,'first');
    Cutoff_0_y(i)=zn1(idx);
    Cutoff_0_z(i)=PercDiff2(i,idx);
end

ciccio_y=movmean(Cutoff_0_y,50);

%h_0=plot3(Cutoff_0_x,Cutoff_0_y,Cutoff_0_z,'-','linewidth',0.1);
%h_0=plot3(Cutoff_0_x,Cutoff_0_y,ciccio_z+0.2,'-','linewidth',0.1);
h_0=plot3(Cutoff_0_x,ciccio_y,zeros(length(za1))+0.01,'k','linewidth',1.5);

text(1.5,2.5,0.3,'Non-Ag Workers: \Delta U^i>0');
text(1.25,0.5,0,'Ag Workers: \Delta U^i<0');
text(0,0.75,0.3,'Non-Ag Workers:');
text(0,0.5,0.2,'\Delta U^i<0');

uistack(h_0,'top');
figure_name = '2D_Change_In_Utility_If_Completely_Remove_Subsidy_GE';
r=150;set(gcf,'PaperPosition',[0 0 1000 800]/r);
% print('-depsc2',figure_name);


%% Counterfactual IV: Remove PSE and add no subsidy (ADD NOTHING), but holding price p fixed at the level with PSE (NP=no subsidy, partial equilibrium), look at welfare of everyone (za,zn)
% Argue that PSE is NOT a subsidy to farmers, but a subsidy to the poor
[NP.yaT,NP.ynT,NP.AgEmpT,NP.NAgEmpT,NP.SubAgT]=deal(nan(T,1));
NP.AgriIndV2=nan(N*N,T);
NP.IncT3=nan(N,N,T); % bring home income
NP.pT=B.pT; % initial guess

for t=1:T
    A=AT(t);
    
    yn=@(p)A*sum((p.*zaV<znV).*znV.*gV);
    ya=@(p)A*sum((p.*zaV>=znV).*zaV.*gV);
    p=NP.pT(t);
    
    NP.AgriIndV=(p*zaV>=znV); 
    NP.AgriIndV2(:,t)=NP.AgriIndV;
    NP.yaT(t)=A*sum(NP.AgriIndV.*zaV.*gV); % in units of crops
    NP.ynT(t)=A*sum((1-NP.AgriIndV).*znV.*gV);      % in units of cars
    
    NP.AgEmpT(t)=sum(NP.AgriIndV.*gV);
    NP.NAgEmpT(t)=sum((1-NP.AgriIndV).*gV);
    NP.IncT3(:,:,t)=A*max(zn2,za2.*p); % income
end
NP.UT3=u2(repmat(permute(B.pT,[2 3 1]),[N N 1]),NP.IncT3); % utility of individual (za,zn) in period t
NP.AvgAgIncT=NP.pT.*NP.yaT./NP.AgEmpT; % in units of cars
NP.AvgNAgIncT=NP.ynT./NP.NAgEmpT; % in units of cars
NP.GDPT=NP.ynT+NP.pT.*NP.yaT;
NP.AgGDPT=NP.pT.*NP.yaT; % Agri GDP

% change in utility of everyone (over the entire time frame) as a function of (za,zn)
figure(51);clf;Period_of_Interest=1;
Utemp2=B.UT3(:,:,Period_of_Interest);
Utemp2(imag(Utemp2)~=0)=nan;
NP.Utemp2=NP.UT3(:,:,Period_of_Interest);
NP.Utemp2(imag(NP.Utemp2)~=0)=nan;
UtilityDiff2=NP.Utemp2(D1,D2)-Utemp2(D1,D2);
surf(za2(D1,D2),zn2(D1,D2),UtilityDiff2); % first period
xlabel z^a;ylabel z^n;zlabel Utility;
xlim([0 za1(N/Frac)]);ylim([0 zn1(N/Frac)]);
figure_name = 'Change_In_Utility_If_Completely_Remove_Subsidy_PE';
r=150;set(gcf,'PaperPosition',[0 0 600 500]/r);
% print('-depsc2',figure_name);


% Agri Emp Share (before and after subsidy removal)
figure(53);clf;
plot(Y,B.AgEmpT*100,'b-',Y,NO.AgEmpT*100,'r--',Y,NP.AgEmpT*100,'magenta-.');grid on;box on;
% legend('PSE','PSE Removed (PE+GE)','PSE Removed, Price Fixed (PE)','location','SE');box on;grid on;
text(2002.6,1.72,'with s_t','fontsize',14);
text(2002.6,1.37,'w/o s_t','fontsize',14);
text(2002.6,1.315,'General Equilibrium','fontsize',13);
text(2002.6,0.57,'w/o s_t','fontsize',14);
text(2002.6,0.5,'Partial Equilibrium','fontsize',13);
xlabel Year;ylabel ('Agriculture Employment Share');
ylim([0.2 1.9]);xlim([CaliStartY CaliEndY]);
figure_name = 'Agri_Emp_Share_NO_NP';
r=150;set(gcf,'PaperPosition',[0 0 600 500]/r);
% print('-depsc2',figure_name);



% change in utility of everyone, as a function of (za,zn), projected into 2D
figure(54);clf;Period_of_Interest=1;
Utemp2=B.UT3(:,:,Period_of_Interest);
Utemp2(imag(Utemp2)~=0)=nan;
NP.Utemp2=NP.UT3(:,:,Period_of_Interest);
NP.Utemp2(imag(NP.Utemp2)~=0)=nan;
UtilityDiff2=NP.Utemp2-Utemp2;
PercUtilityDiff2=UtilityDiff2./(-(Utemp2-max(Utemp2(:))))*100;
% h=surf(za2(D1,D2),zn2(D1,D2),PercUtilityDiff2(D1,D2)); % first period
h=surf(za2,zn2,PercUtilityDiff2); % first period
set(h,'edgecolor','none');
colorbar;
caxis([-10 0]);
xlabel z^a;ylabel z^n;zlabel 'Percentage Change in Utility';title('Partial Equilibrium Effect');
% view([0 90]);
view([-45 40]);
xlim([0 za1(N/4)]);ylim([0 za1(N/4)]);zlim([-8 0.5]);
Period_of_Interest=1;
Cutofftemp2=reshape(NP.AgriIndV2(:,Period_of_Interest),[N N]); % agriculture occupation indicator
Cutoff_x=za1;
Cutoff_y=nan(N,1);
Cutoff_z=nan(N,1);
for i=1:N % over zn
   idx=find(Cutofftemp2(i,:)==1,1,'last'); % the highest zn such that households still stay in agriculture
   Cutoff_y(i)=zn1(idx);
   Cutoff_z(i)=PercUtilityDiff2(i,idx);
end
ciccio_y=movmean(Cutoff_y,20);
ciccio_z=movmean(Cutoff_z,50);
hold on;
h=plot3(Cutoff_x,ciccio_y,Cutoff_z,'r','linewidth',1.5);
uistack(h,'top');
Cutoff_0_x=za1;
Cutoff_0_y=nan(N,1);
Cutoff_0_z=nan(N,1);
for i=1:N
    temp=UtilityDiff2(i,:);
    idx=find(temp<0,1,'last');
    Cutoff_0_y(i)=zn1(idx);
    Cutoff_0_z(i)=PercUtilityDiff2(i,idx);
end
clear ciccio_y ciccio_z
ciccio_y=movmean(Cutoff_0_y,20);
ciccio_z=movmean(Cutoff_0_z,50);

h_0=plot3(Cutoff_0_x,ciccio_y,zeros(length(za1)),'k','linewidth',1.5);
uistack(h_0,'top');

text(1.5,2.5,0.3,'Non-Ag Workers: \Delta U^i>0');
text(1.25,0.1,-4,'Ag Workers: \Delta U^i<0');
figure_name = '2D_Change_In_Utility_If_Completely_Remove_Subsidy_PE';
r=150;set(gcf,'PaperPosition',[0 0 1000 800]/r);
% print('-depsc2',figure_name);


%% Counterfactual VII: Trump Tariff Short Run (TS: Trump Short)
% Modeling China's retaliation tariff as a reduction in gam0 (international demand)
% In the SR, ag and non-ag career choice is fixed at the Benchmark level
% Ag price tanks, calculate the drop in price and compensation for farmers.

TS.gam0=gam0*(1+25/100)^(-gam1); % reduce international demand due to 25% Chinese tariff on US crops

% Calculate the time series of consumption, employment, GDP from starting year to ending year
[TS.yaT,TS.ynT]=deal(nan(T,1));
TS.ImpliedIncT=nan(N,N,T); % implied income to bring utility back to Benchmark level
TS.pT=B.pT; % initial guess
pLB=1e-3;
pUB=pBase*1.5;
TS.Total_Farmer_CompensationT=nan(T,1); % total lump-sum Trump compensation to farmers due to tariff (as a percentage of Agri GDP)


for t=1:T
    s=sT(t);
    A=AT(t);
    
    if Y(t)<=CaliEndY
        g0=gam0;
    else 
        g0=TS.gam0;
    end
    
    yn=B.ynT(t); % in SR, output remains at Benchmark level
    ya=B.yaT(t); % in SR, output remains at Benchmark level
    ca=@(p)ya.*(1-g0*p.^(-gam1)); % domestic consumption of agri product
    
    p=fzero(@(p)yn-nu*p*(ca(p)-abar),[pLB pUB],optimset('display','iter','tolX',pLR/1e3)); % p solves FOC of households
    pLB=max(0,p*0.5); % new lower bound
    pUB=p*1.3; % new upper bound
    
    TS.AgriIndV=B.AgriIndV;                   % indicator for working in Ag: use Benchmark level indicator
    TS.yaT(t)=ya;                                       % agri output in units of crops
    TS.ynT(t)=yn;                                       % non-agri output in units of cars
    TS.pT(t)=p;                                           % relative price of Ag goods in year t
    
    TS.AgEmpT(t)=sum(TS.AgriIndV.*gV);                  % total/proportion of employment in agri in period t
    TS.NAgEmpT(t)=sum((1-TS.AgriIndV).*gV);         % total/proportion of employment in non-agri in period t 
    
    TS.ImpliedIncT3(:,:,t)=p*abar+exp((B.UT3(:,:,t)+log(p)+(1+nu)*log(1+nu)-nu*log(nu))/(1+nu));  % implied bring-home income needed to generate utility level of Benchmark model
end

% solve for compensation to farmers required to bring their utility back to Benchmark level
tau_diff=1000;
tau_iter=1;
tau_iterMax=1e3;

tauTGuess=B.tauT; % use benchmark level tax rate as initial guess
tauTUpdate=nan(T,1);

while tau_diff>1e-7 && tau_iter<tau_iterMax
    tau_iter=tau_iter+1;
    
    for t=1:T
        AgriInd_Benchmark=(B.pT(t)*(1+sT(t))*za2>=zn2); % indicator for a worker being in agri sector before Trump's tariff
        
        Compensation2=TS.ImpliedIncT3(:,:,t).*AgriInd_Benchmark...
            -AT(t)*(1-tauTGuess(t))*TS.pT(t)*(1+sT(t))*AgriInd_Benchmark.*za2; % only compensation whoever used to be farmer before tariff
        
        tauTUpdate(t)=(sT(t)*TS.pT(t)*TS.yaT(t)+sum(Compensation2(:).*gV))/(TS.yaT(t)*TS.pT(t)+TS.ynT(t));
        
        TS.Farmer_CompensationT(:,:,t)=Compensation2;
        TS.Total_Farmer_CompensationT(t)=sum(sum(Compensation2(:).*gV))/(TS.pT(t)*TS.yaT(t))*100; % total compensation as fraction of (new) agriculture GDP
    end
    tau_diff=max(abs(tauTGuess-tauTUpdate))
    tauTGuess=(tauTGuess+tauTUpdate)/2;
end
TS.tauT=tauTGuess;  % proportional tax rate on everyone to finance PSE spending and compensation to farmers

% Relative Price of Agri Output (before and after subsidy removal)
figure(81);clf;
plot(Y,B.pT,'b-',Y,TS.pT,'r--');
legend('Benchmark (PSE): p','Trump Tariff (SR): p','location','SW');box on;grid on;
xlabel Year;ylabel ('Relative Price of Agriculture Output');
figure_name = 'Relative_Price_Agri_Output_Trump_SR';
r=150;set(gcf,'PaperPosition',[0 0 600 500]/r);
% print('-depsc2',figure_name);

% total compensation to ex-farmer required in the SR
figure(83);clf;
plot(Y,TS.Total_Farmer_CompensationT,'bo-');
box on;grid on;
xlabel Year;ylabel ('Percentage of Agriculture GDP');
figure_name = 'Total_Compensation_to_Farmers_Trump_Tariff_SR';
r=150;set(gcf,'PaperPosition',[0 0 600 500]/r);
% print('-depsc2',figure_name);



%% Counterfactual VI: Trump Tariff Long Run (TL: Trump Long) 
% Modeling China's retaliation tariff as a reduction in gam0 (international demand)
% In the SR, ag and non-ag career choice is fixed at the Benchmark level. Ag price tanks, calculate the drop in price and compensation for farmers.
% In the LR, ag and non-ag career choice is flexible. Farmers relocate to non-agri sector, resulting in a reduction in supply of ag goods. Price of agri doesn't drop as much and compensation required is smaller. 
[TL.yaT,TL.ynT,TL.AgEmpT,TL.NAgEmpT]=deal(nan(T,1));
TL.Farmer_CompensationT=nan(N,N,T); % lump-sum Trump compensation to farmers due to tariff
TL.Total_Farmer_CompensationT=nan(T,1); % total lump-sum Trump compensation to farmers due to tariff (as a percentage of Agri GDP)
TL.ImpliedIncT3=nan(N,N,T); % bring-home income required to generate the utility level of Benchmark
TL.pT=nan(T,1); % LR price after tariff war
pLB=1e-3;
pUB=pBase*1.5;

for t=1:T
    s=sT(t);
    A=AT(t);
    
    if Y(t)<=CaliEndY
        g0=gam0;
    else 
        g0=TS.gam0;
    end

    yn=@(p)A*sum((p.*(1+s).*zaV<znV).*znV.*gV); % non agri output
    ya=@(p)A*sum((p.*(1+s).*zaV>znV).*zaV.*gV); % agri output
    ca=@(p)ya(p).*(1-g0*p.^(-gam1)); % domestic consumption of agri product
    
    p=fzero(@(p)yn(p)-nu*p*(ca(p)-abar),[pLB pUB],optimset('display','iter','tolX',pLR/1e3)); % p solves FOC of households
    pLB=max(0,p*0.5); % new lower bound
    pUB=p*1.3; % new upper bound
    
    TL.pT(t)=p;
    TL.AgriIndV=(p*(1+s)*zaV>=znV);                   % indicator for working in Ag
    TL.yaT(t)=A*sum(TL.AgriIndV.*zaV.*gV);             % agri output in units of crops
    TL.ynT(t)=A*sum((1-TL.AgriIndV).*znV.*gV);      % non-agri output in units of cars

    TL.AgEmpT(t)=sum(TL.AgriIndV.*gV);                  % total/proportion of employment in agri in period t
    TL.NAgEmpT(t)=sum((1-TL.AgriIndV).*gV);         % total/proportion of employment in non-agri in period t
    
    TL.ImpliedIncT3(:,:,t)=p*abar+exp((B.UT3(:,:,t)+log(p)+(1+nu)*log(1+nu)-nu*log(nu))/(1+nu));  % bring-home income needed to generate utility level as Benchmark level
end

tau_diff=1000;
tau_iter=1;
tau_iterMax=1e3;

tauTGuess=B.tauT; % use benchmark level tax rate as initial guess
tauTUpdate=nan(T,1);

while tau_diff>1e-7 && tau_iter<tau_iterMax
    tau_iter=tau_iter+1;
    
    for t=1:T
        AgriInd_Benchmark=(B.pT(t)*(1+sT(t))*za2>=zn2); % indicator for a worker being in agri sector before Trump's tariff
        Compensation2=TL.ImpliedIncT3(:,:,t).*AgriInd_Benchmark...
            -AT(t)*(1-tauTGuess(t))*TL.pT(t)*(1+sT(t))*AgriInd_Benchmark.*za2; % only compensation whoever used to be farmer before tariff
        tauTUpdate(t)=(sT(t)*TL.pT(t)*TL.yaT(t)+sum(Compensation2(:).*gV))/(TL.yaT(t)*TL.pT(t)+TL.ynT(t));
        TL.Farmer_CompensationT(:,:,t)=Compensation2;
        TL.Total_Farmer_CompensationT(t)=sum(sum(Compensation2(:).*gV))/(TL.pT(t)*TL.yaT(t))*100; % total compensation as fraction of (new) agriculture GDP
    end
    tau_diff=max(abs(tauTGuess-tauTUpdate))
    tauTGuess=(tauTGuess+tauTUpdate)/2;
end
TL.tauT=tauTGuess;  % proportional tax rate on everyone to finance PSE spending and compensation to farmers



% Relative Price of Agri Output (before and after subsidy removal)
figure(91);clf;
plot(Y,B.pT,'b-',Y,TS.pT,'m-.',Y,TL.pT,'r--');
legend('No Tariff','Tariff: Short Run','Tariff: Long Run','location','SW');box on;grid on;
xlabel Year;ylabel ('Relative Price of Agriculture Output');
figure_name = 'Relative_Price_Agri_Output_Trump_SR_LR';
r=150;set(gcf,'PaperPosition',[0 0 600 500]/r);
% print('-depsc2',figure_name);

for i=1:17
    TS.Total_Farmer_CompensationT2(i)=0;
end

for i=18:20
    TS.Total_Farmer_CompensationT2(i)=TS.Total_Farmer_CompensationT(i);
end

% total compensation to ex-farmer required in the SR and LR (as % of agri GDP)
figure(93);clf;
plot(Y(13:20),TS.Total_Farmer_CompensationT2(13:20),'m-.',Y(13:20),TL.Total_Farmer_CompensationT(13:20),'r--');
box on;grid on;
legend('Short Run','Long Run','location','NW');
xlabel Year;ylabel ('% of Agriculture GDP');
figure_name = 'Total_Compensation_to_Farmers_Trump_Tariff_SR_LR';
r=150;set(gcf,'PaperPosition',[0 0 600 500]/r);
% print('-depsc2',figure_name);


for i=1:17
    TS.pT2(i)=B.pT(i);
end

for i=18:20
    TS.pT2(i)=TS.pT(i);
end

figure(99);clf;
plot(Y(13:20),B.pT(13:20),'b-',Y(13:20),TS.pT2(13:20),'magenta-.');
box on;grid on;
hold
plot(Y(13:20),TL.pT(13:20),'r--');
xlabel Year;
ylabel ('p^a')
legend('No Tariff','Short Run','Long Run','location','NW');
figure_name = ('patariffs');
r=150;set(gcf,'PaperPosition',[0 0 600 500]/r);
% print('-depsc2',figure_name);


figure(991);clf;
plot(Y,B.pT,'b-',Y,TS.pT2,'magenta-.');
box on;grid on;
hold
plot(Y,TL.pT,'r--');
xlim([2006 2020])
xlabel Year;
ylabel ('p^a')
legend('No Tariff','Short Run','Long Run','location','SW');
figure_name = ('patariffs');
r=150;set(gcf,'PaperPosition',[0 0 600 500]/r);
% print('-depsc2',figure_name);


%% Statistics used for Introduction of Draft

% statistics for draft
Yr=2006; % baseline year (2006 is selected because it is right before Financial Crisis)
YB=(Y==Yr); % baseline year indicator

fprintf('\n')
fprintf('\n')
fprintf('\n')
disp(strcat('Baseline year PSE=',num2str(sT(YB)*100),'%'));fprintf('\n');
disp(strcat('Price change without PSE in base year...',num2str(((LS.pT(YB)-B.pT(YB))./B.pT(YB)) *100),'%'));fprintf('\n');
disp(strcat('Average Price change without PSE...',num2str(mean((LS.pT-B.pT)./B.pT) *100),'%'));fprintf('\n');

disp(strcat('Wage loss after GE effect...',num2str(-(LS.pT(YB)-B.pT(YB)*(1+sT(YB)))/(B.pT(YB)*(1+sT(YB)))*100),'%'));fprintf('\n');
disp(strcat('Replacing PSE with Lump-Sum Transfer leads to  ',num2str((LS.pT(YB)-B.pT(YB))/B.pT(YB)*100),'% increase in Agri Price'));fprintf('\n');
disp(strcat('Avg Ag worker income drops by ',num2str((B.AvgAgInc_Include_Subsidy_T(YB)-LS.AvgAgIncT(YB))/B.AvgAgInc_Include_Subsidy_T(YB)*100),'%'));fprintf('\n');

Benchmark_APG=B.AvgAgInc_Include_Subsidy_T./B.AvgNAgIncT;
Counterfactual_APG=NO.AvgAgIncT./NO.AvgNAgIncT;
disp(strcat('APG Data (base year)... ',num2str(Benchmark_APG(YB))));fprintf('\n');
disp(strcat('APG No ta (base year)... ',num2str((Counterfactual_APG(YB)))));fprintf('\n');

disp(strcat('GE Fiscal Saving from Replacing PSE with Lump-Sum Transfer =... ',num2str(LS.RevT(YB)./LS.AgGDPT(YB)*100),'%... of Ag GDP'));fprintf('\n');
disp(strcat('PE Fiscal Saving from Replacing PSE with Lump-Sum Transfer =... ',num2str(PE.RevT(YB)./PE.AgGDPT(YB)*100),'%... of Ag GDP'));fprintf('\n');

disp(strcat('GE Average Fiscal Saving from Replacing PSE with Lump-Sum Transfer =... ',num2str(mean(LS.RevT./LS.AgGDPT)*100),'%... of Ag GDP'));fprintf('\n');
disp(strcat('PE Average Fiscal Saving from Replacing PSE with Lump-Sum Transfer =... ',num2str(mean(PE.RevT./PE.AgGDPT)*100),'%... of Ag GDP'));fprintf('\n');


Yr=2018;
YB=(Y==Yr); % baseline year indicator

disp(strcat('China Retalitatory tariffs lead to...',num2str((B.pT(YB)-TS.pT(YB))/B.pT(YB)*100),'%.. drop in relative price of ag goods in the Short Run')); fprintf('\n');
disp(strcat('China Retalitatory tariffs lead to... ',num2str((B.pT(YB)-TL.pT(YB))/B.pT(YB)*100),'%.. drop in relative price of ag goods in the Long Run')); fprintf('\n');

disp(strcat('China Retalitatory tariffs require... ',num2str(TS.Total_Farmer_CompensationT(YB)),'%.. of Ag VA as compensation in the Short Run')); fprintf('\n');
disp(strcat('China Retalitatory tariffs require... ',num2str(TL.Total_Farmer_CompensationT(YB)),'%.. of Ag VA as compensation in the Long Run')); fprintf('\n');

figure(999);clf;
plot(Y,Benchmark_APG,'b-',Y,Counterfactual_APG,'r--');
box on;grid on;
legend('Benchmark (with s_t)','Counterfactual (w/o s_t)','location','SE');
xlabel Year;ylabel ('APG');
figure_name = 'APG';
r=150;set(gcf,'PaperPosition',[0 0 600 500]/r);
xlim([2001 2017])
% print('-depsc2',figure_name);
