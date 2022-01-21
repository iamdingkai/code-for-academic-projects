function [GapAvgWage,pBase,za1,zn1,za2,zn2,g2,zaV,znV,gV] ...
    = Fun_Cali_rho(rho,AgEmpBase,sBase,AvgWageRatioBase,thetaa,thetan,N)

global pBase_Guess;

% discretizing the distribution of (za,zn)
C=@(u,v)-1./rho.*log(1+((exp(-rho*u)-1).*(exp(-rho*v)-1))./(exp(-rho)-1)); % joint distribution of productivity draws
F=@(za)exp(-za.^(-thetaa));
H=@(zn)exp(-zn.^(-thetan));
G=@(za,zn)C(F(za),H(zn));               % CDF
f=@(za)thetaa*F(za).*za.^(-1-thetaa);
h=@(zn)thetan*H(zn).*zn.^(-1-thetan);
C_uv=@(u,v)-rho*exp(-rho*(u+v)).*(exp(-rho)-1)./((exp(-rho)-1)+(exp(-rho*u)-1).*(exp(-rho*v)-1)).^2;
g=@(za,zn)C_uv(F(za),H(zn)).*f(za).*h(zn); % PDF
% calculate mean and standard deviation of za and zn
sp=30;
mu_za=gamma(1-1/thetaa);
mu_zn=gamma(1-1/thetan);
sd_za=sqrt(gamma(1-2/thetaa)-(gamma(1-1/thetaa))^2);
sd_zn=sqrt(gamma(1-2/thetan)-(gamma(1-1/thetan))^2);
% grid for za and zn
za1=linspace(1e-10,mu_za+sp*sd_za,N+2)';
za1=za1(3:N+2);
zn1=linspace(1e-10,mu_zn+sp*sd_zn,N)';
% joined grids for za and zn
[za2,zn2]=ndgrid(za1,zn1); 
g2=g(za2,zn2);
g2=g2./sum(g2(:));
% reshape joined grids into a vector for faster calculation
zaV=za2(:);
znV=zn2(:);
gV=g2(:);

%%  Base parameteres
AgEmpBase_handle=@(p)sum(gV.*(p*(1+sBase).*zaV>=znV)); % percentage emp in agri
o=optimset('Display','on');
pBase=fzero(@(p)100*(AgEmpBase_handle(p)-AgEmpBase),pBase_Guess,o);
pBase_Guess=pBase;

AgriIndV=(pBase*(1+sBase).*zaV>=znV);

AvgWageAg=pBase*(1+sBase)*sum(gV.*zaV.*AgriIndV)/AgEmpBase;
AvgWageNAg=sum(gV.*znV.*(1-AgriIndV))/(1-AgEmpBase);

AvgWageRatioModel=AvgWageAg/AvgWageNAg;

GapAvgWage=AvgWageRatioBase-AvgWageRatioModel;

end

