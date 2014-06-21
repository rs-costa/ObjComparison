function dist=eclDistance(glob,x)
noOfRatios=length(glob.exp_flux);
W=zeros(noOfRatios,noOfRatios);
for i=1:noOfRatios
 	W(i,i)=(1/sum(1./glob.var_exp))*(1/glob.var_exp(i));
end
dist=sqrt((x(glob.indices)-glob.exp_flux)'*W*(x(glob.indices)-glob.exp_flux));
end