args <- commandArgs(trailingOnly = TRUE)

# Bakteri katvuste tabel
bac_file <- args[1]
plasmid_file <- args[2]
read_length <- as.numeric(args[3])
word <- as.numeric(args[4])
coverage_variation <- as.numeric(args[5])
output_name <- args[6]

use_legacy = 0

# kasutusn2ited on toodud peale funktsiooni defineerimist - andmete paiknemiskataloog tuleb muuta!!!


# Funktsiooni definitsiooni algus *********************************************************************************************

testiPlasmiidi2=function(k, k1, n1, k2, n2, readi_pikkus, pr=0.05){
# Argumendid
# k - kmeri pikkus
# k1 - bakteri katvused
# n1 - mitu bakteri k-meeri oli sellise katvusega
# k2 - plasmiidi katvused 
# n2 - mitu plasmiidi k-meeri oli sellise katvusega
# readi_pikkus - kui pikki lugemeid on kasutatud
# pr - kui suur on lubatud sequencing bias (mis tekib n2iteks cg-osakaalu erinevusest)


par_to_p = function(x){
 # 22rmised t6en2osused 0...0.3  keskmine 1..0.4
 p00=0.3*exp(x)/(1+exp(x))
 p=c(p00[1], 1-p00[1]-p00[2], p00[2])
 return(p)
}

p_to_par = function(x){
  x[x==0]=1e-7
  xx=log( 1/0.3* x[c(1,3)] /(1-  (1/0.3) *x[c(1,3)]) )
  return( xx )
}

par_to_kesk=function(x, keskm0=keskm0){
# 0.6*keskm0 ... 1.4*keskm0
# return(0.8*keskm0*exp(x)/(1+exp(x))+0.6*keskm0)
return(exp(x))
}



kesk_to_par=function(x, keskm0=keskm0){
# 0.6*keskm0 ... 1.4*keskm0
# abi=(x-0.6*keskm0)/(0.8*keskm0)
# return( log(abi/(1-abi) ) )
return(log(x))
}

# p_to_par=c(0.2,0.4, 0.2)
# Keskmist t6en2osust ignoreeritakse...
# par_to_p(  p_to_par(c(0.2,0.4, 0.24))  )

ll=function(x, k, n, keskm0=NA, yle=NULL, ypiir){
# print ("Funktsioon ll")
#  print (x)
  # pi - 0x 1x 2x 3x esinevate k-meeride osakaalud
  # x=c(1,1, 12, 5); k=a5$katvus; n=a5$n
  # x=c(p_to_par(c(0.05, 0.9, 0.05)), kesk_to_par(12), 5); k=a1$katvus; n=a1$n; yle=20; keskm0=13
  # x=c(p_to_par(c(0.05, 0.9, 0.05)), kesk_to_par(12), 5); k=a5$katvus; n=a5$n; yle=20; keskm0=13
# tul2=optim(c(p_to_par(c(0.05, 0.9, 0.05)), kesk_to_par(keskm00_plasmiid, keskm00) ), ll, k=k2, n=n2, control=list(fnscale=-1), keskm0=keskm00, yle=exp(tul$par[4]) )
  # x=c(p_to_par(c(0.05, 0.9, 0.05)), kesk_to_par(keskm00_plasmiid, keskm00)); k=a5$katvus; n=a5$n; yle=exp(tul$par[4]); keskm0=keskm00



  pi=par_to_p( x[1:2] )
  
  # keskmine katvus
  keskm=par_to_kesk(x[3], keskm0)

  # ylehajuvusparameeter
  if (is.null(yle)) s=exp(x[4]) else s=yle

  if (keskm<=0) return(-Inf) else {
  # yhe vaatluse kontributsioon
#  l1=log(dpois(k, lambda=0.1)*pi[1]+dpois(k, lambda=keskm)*pi[2]+dpois(k, lambda=keskm*2)*pi[3]+dpois(k, lambda=keskm*3)*pi[4])

# Arvutuste t2psuse t6stmine
# ln(a + b) = ln{exp[ln(a) - ln(b)] + 1} + ln(b)

# 1 v6i 0 t6en2osus (Roosaare failides 0-d sageli puudu)
l1_1 =pnbinom(1, mu=0.1, size=s, log=T)+log(pi[1])
l1_2 =pnbinom(1, mu=keskm, size=s, log=T)+log(pi[2])
l1_3 =pnbinom(1, mu=keskm*2, size=s, log=T)+log(pi[3])
# l1_4 =pnbinom(1, mu=keskm*3, size=s, log=T)+log(pi[4])

# J2rjekord selline, sest exp(l1_1)>0 ka siis, kui keskm~l6pmatus
l1_a = log(exp(l1_2-l1_1)+1)+l1_1
# l1_b = log(exp(l1_3-l1_4)+1)+l1_4
l1_b = l1_3
l1_1v = log(exp(l1_b-l1_a)+1)+l1_a

# V22rtusest 2*keskm0 suurema v22rtuse n2gemise t6en2osus
l1_1 =pnbinom( ypiir, mu=0.1, size=s, log=T, lower.tail=FALSE)+log(pi[1])
l1_2 =pnbinom(ypiir, mu=keskm, size=s, log=T, lower.tail=FALSE)+log(pi[2])
l1_3 =pnbinom(ypiir, mu=keskm*2, size=s, log=T, lower.tail=FALSE)+log(pi[3])
# l1_4 =pnbinom(ypiir, mu=keskm*3, size=s, log=T, lower.tail=FALSE)+log(pi[4])
l1_a = log(exp(l1_1-l1_2)+1)+l1_2
# l1_b = log(exp(l1_3-l1_4)+1)+l1_4
l1_b = l1_3
l1_1s = log(exp(l1_a-l1_b)+1)+l1_b


# vahepealsed
ind= 1<k & k<ypiir
k0=k[ind]

l1_1 =dnbinom(k0, mu=0.1, size=s, log=T)+log(pi[1])
l1_2 =dnbinom(k0, mu=keskm, size=s, log=T)+log(pi[2])
l1_3 =dnbinom(k0, mu=keskm*2, size=s, log=T)+log(pi[3])
# l1_4 =dnbinom(k0, mu=keskm*3, size=s, log=T)+log(pi[4])
l1_a = log(exp(l1_1-l1_2)+1)+l1_2
# l1_b = log(exp(l1_3-l1_4)+1)+l1_4
l1_b = l1_3
l1 = log(exp(l1_a-l1_b)+1)+l1_b

  # Valimi log-t6ep2ra
# tricube f(u)=70/81*(1-|u|^3)^3 range |u|<=1
#  lisa=0
#  if (tricube) { u=abs(4*(keskm-keskm0)/keskm0); lisa=-Inf; if(u<1) lisa=log(70/81*(1-u^3 )^3 )  }
  l2=sum(l1*n[ind]) + l1_1v*sum(n[k<=1])+l1_1s*sum(n[k>ypiir])
  if (is.na(l2) | abs(sum(pi)-1)>1e-14 ) l2=-Inf
  return(l2)
}
}



ll_var2 = function(x, k, n, keskmine, keskm0=NA, yle=NULL, protsent, ypiir){
# protsent - n2itab, kui mitme protsendi v6rra v6ib keskmine katvus erineda

  # pi - 0x 1x 2x 3x esinevate k-meeride osakaalud
  # x=c(1,1,1,12); k=a5$katvus; n=a5$n
  pi=par_to_p( x[1:2] )


  # keskmine katvus
  # lubame k6rvalekaldeid keskm0-st vahemikus    keskm0*(1-protsent) .... keskm0*(1+protsent)
  if (is.na(x[3])) keskm=keskmine else keskm= 2*protsent*keskmine*exp(x[3])/(1+exp(x[3]))+ (1-protsent)*keskmine


  # ylehajuvusparameeter
  if (is.null(yle)) s=exp(x[4]) else s=yle


  # yhe vaatluse kontributsioon

# Arvutuste t2psuse t6stmine
# ln(a + b) = ln{exp[ln(a) - ln(b)] + 1} + ln(b)

# 1 v6i 0 t6en2osus (Roosaare failides 0-d sageli puudu)
l1_1 =pnbinom(1, mu=0.1, size=s, log=T)+log(pi[1])
l1_2 =pnbinom(1, mu=keskm, size=s, log=T)+log(pi[2])
l1_3 =pnbinom(1, mu=keskm*2, size=s, log=T)+log(pi[3])
# l1_4 =pnbinom(1, mu=keskm*3, size=s, log=T)+log(pi[4])
l1_a = log(exp(l1_1-l1_2)+1)+l1_2
# l1_b = log(exp(l1_3-l1_4)+1)+l1_4
l1_b = l1_3
l1_1v = log(exp(l1_a-l1_b)+1)+l1_b

# V22rtusest 2*keskm0 suurema v22rtuse n2gemise t6en2osus
l1_1 = pnbinom( ypiir, mu=0.1, size=s, log=T, lower.tail=FALSE)+log(pi[1])
l1_2 = pnbinom(ypiir, mu=keskm, size=s, log=T, lower.tail=FALSE)+log(pi[2])
l1_3 = pnbinom(ypiir, mu=keskm*2, size=s, log=T, lower.tail=FALSE)+log(pi[3])
# l1_4 = pnbinom(ypiir, mu=keskm*3, size=s, log=T, lower.tail=FALSE)+log(pi[4])
l1_a = log(exp(l1_1-l1_2)+1)+l1_2
# l1_b = log(exp(l1_3-l1_4)+1)+l1_4
l1_b = l1_3
l1_1s = log(exp(l1_a-l1_b)+1)+l1_b


# vahepealsed
ind= 1<k & k<ypiir
k0=k[ind]

l1_1 =dnbinom(k0, mu=0.1, size=s, log=T)+log(pi[1])
l1_2 =dnbinom(k0, mu=keskm, size=s, log=T)+log(pi[2])
l1_3 =dnbinom(k0, mu=keskm*2, size=s, log=T)+log(pi[3])
# l1_4 =dnbinom(k0, mu=keskm*3, size=s, log=T)+log(pi[4])
l1_a = log(exp(l1_1-l1_2)+1)+l1_2
# l1_b = log(exp(l1_3-l1_4)+1)+l1_4
l1_b = l1_3
l1 = log(exp(l1_a-l1_b)+1)+l1_b


  # Valimi log-t6ep2ra
#  lisa=0
#  if (tricube) { u=abs(4*(keskm-keskm0)/keskm0); lisa=-Inf; if(u<1) lisa=log(70/81*(1-u^3 )^3 )  }

  # l2=sum(l1*n)+lisa
  l2=sum(l1*n[ind]) + l1_1v*sum(n[k<=1])+l1_1s*sum(n[k>ypiir])
  if (is.na(l2) | abs(sum(pi)-1)>1e-14 ) l2=-Inf
  return(l2)
}


# bakter
# ...........

keskm00=sum(k1*prop.table(n1))

# ylempiir - millise v22rtuseni vaadeldakse
# ypiir=ceiling(keskm00*2)
ypiir=ceiling(keskm00*5)

tul=optim(c(p_to_par(c(0.05, 0.9, 0.05)),  kesk_to_par(keskm00, keskm00), 4.5 ), ll, k=k1, n=n1, control=list(fnscale=-1), keskm0=keskm00, ypiir=ypiir)
tul=optim( tul$par, ll, k=k1, n=n1, control=list(fnscale=-1), keskm0=keskm00, ypiir=ypiir)
tul=optim( tul$par, ll, k=k1, n=n1, control=list(fnscale=-1), keskm0=keskm00, ypiir=ypiir)
tul=optim( tul$par, ll, k=k1, n=n1, control=list(fnscale=-1), keskm0=keskm00, ypiir=ypiir)
tul=optim( tul$par, ll, k=k1, n=n1, control=list(fnscale=-1), keskm0=keskm00, ypiir=ypiir)
tul=optim( tul$par, ll, k=k1, n=n1, control=list(fnscale=-1), keskm0=keskm00, ypiir=ypiir)
tul=optim( tul$par, ll, k=k1, n=n1, control=list(fnscale=-1), keskm0=keskm00, ypiir=ypiir)

count=0
while (tul$convergence!=0 & count<20){
  tul=optim( tul$par, ll, k=k1, n=n1, control=list(fnscale=-1), keskm0=keskm00, ypiir=ypiir)
  count=count+1
}

 pii1=par_to_p( tul$par[1:2] )

if (pii1[2]<0.65) warning("Wrong bacteria or convergence problem? Estimated proportion of bacterial k-mers which are present and unique is too small! (<0.65)")


# plasmiid
# ...........

keskm00_plasmiid=sum(k2*prop.table(n2))

# tul2=optim(c(2, 1, 1, (sum(katvus*prop.table(nn2)))  ), ll, n=nn2, control=list(fnscale=-1), ypiir=ypiir)
tul2=optim(c(p_to_par(c(0.05, 0.9, 0.05)), kesk_to_par(keskm00_plasmiid, keskm00) ), ll, k=k2, n=n2, control=list(fnscale=-1), keskm0=keskm00, yle=exp(tul$par[4]), ypiir=ypiir )
tul2=optim( tul2$par, ll, k=k2, n=n2, control=list(fnscale=-1), keskm0=keskm00, yle=exp(tul$par[4]), ypiir=ypiir )
tul2=optim( tul2$par, ll, k=k2, n=n2, control=list(fnscale=-1), keskm0=keskm00, yle=exp(tul$par[4]), ypiir=ypiir)
tul2=optim( tul2$par, ll, k=k2, n=n2, control=list(fnscale=-1), keskm0=keskm00, yle=exp(tul$par[4]), ypiir=ypiir)
tul2=optim( tul2$par, ll, k=k2, n=n2, control=list(fnscale=-1), keskm0=keskm00, yle=exp(tul$par[4]), ypiir=ypiir)

count=0
while (tul2$convergence!=0 & count<20){
  tul2=optim( tul2$par, ll, k=k2, n=n2, control=list(fnscale=-1), keskm0=keskm00, yle=exp(tul$par[4]), ypiir=ypiir)
  count=count+1
}

 pii2=par_to_p( tul2$par[1:2] )
# if (pii2[2]<0.8) warning("Vale plasmiid v6i koondumisprobleem????")


if (pr<0.001)  tul3=optim(  c(tul2$par[1:2]), ll_var2, k=k2, n=n2, keskmine=par_to_kesk(tul$par[3], keskm00) ,control=list(fnscale=-1), keskm0=keskm00, yle=exp(tul$par[4]), protsent=pr, ypiir=ypiir) else tul3=optim(  c(tul2$par[1:3]), ll_var2, k=k2, n=n2, keskmine=par_to_kesk(tul$par[3], keskm00) ,control=list(fnscale=-1), keskm0=keskm00, yle=exp(tul$par[4]), protsent=pr, ypiir=ypiir)
tul3=optim(  tul3$par, ll_var2, k=k2, n=n2, keskmine=par_to_kesk(tul$par[3], keskm00) ,control=list(fnscale=-1), keskm0=keskm00, yle=exp(tul$par[4]), protsent=pr, ypiir=ypiir)

count=0
while (tul3$convergence!=0 & count<20){
  tul3=tul3=optim(  tul3$par, ll_var2, k=k2, n=n2, keskmine=par_to_kesk(tul$par[3], keskm00) ,control=list(fnscale=-1), keskm0=keskm00, yle=exp(tul$par[4]), protsent=pr, ypiir=ypiir)
  count=count+1
}

 pii3=par_to_p( tul3$par[1:2] )

  ajutkeskm=par_to_kesk(tul$par[3], keskm00)
  if (is.na(tul3$par[3])) keskm_integreeritud =  ajutkeskm else keskm_integreeritud = 2*pr*ajutkeskm*exp(tul3$par[3])/(1+exp(tul3$par[3]))+ (1-pr)*ajutkeskm


like3a=ll(tul2$par, k=k2, n=n2, yle=exp(tul$par[4]), keskm0=keskm00, ypiir=ypiir)
like3b=ll_var2(tul3$par, k=k2, n=n2, keskmine=par_to_kesk(tul$par[3], keskm00), keskm0=keskm00, yle=exp(tul$par[4]), protsent=pr, ypiir=ypiir)

# readi pikkus 100

teststat2 = 2*(like3a-like3b)/(readi_pikkus-k+1)
if (use_legacy){
  teststat2 = 2*(like3a-like3b)/(readi_pikkus)  
}

pvalue2 = 1-pchisq(teststat2, 1)






# vaba Plasmiid - arvutused katvuste erinevuse hindamiseks (pole vajalik testimiseks)



# ylempiir - millise v22rtuseni vaadeldakse
ypiir_plasmiid=ceiling(keskm00_plasmiid*5)

tul4=optim(c(p_to_par(c(0.05, 0.9, 0.05)),  kesk_to_par(keskm00_plasmiid, keskm00_plasmiid), 4.5 ), ll, k=k2, n=n2, control=list(fnscale=-1), keskm0=keskm00_plasmiid, ypiir=ypiir_plasmiid)

count=0
while (tul4$convergence!=0 & count<20){
  tul4=optim( tul4$par, ll, k=k2, n=n2, control=list(fnscale=-1), keskm0=keskm00_plasmiid, ypiir=ypiir_plasmiid)
  count=count+1
}

 pii4=par_to_p( tul4$par[1:2] )

if (pii4[2]<0.65) warning("Plasmid coverage estimation failed???? Missing or wrong plasmid? Estimated proportion of plasmid k-mers which are present and unique is too small (<0.65)!")





koondus=!(tul$convergence!=0 | tul2$convergence!=0 | tul3$convergence!=0 | pii1[2]<0.8)

if (!koondus) { 
    warning("Convergence problem! Estimated parameters and test results are unreliable!")
      #warning(paste("Tehniline probleemi kirjeldus:  koondumine 1:", tul$convergence , 
    #" koondumine 2:", tul2$convergence, " koondumine 3:", tul3$convergence, "bakteri normaalsete k-meeride osakaal:", pii1[2]))
  }

tulem=list(teststat=teststat2, pvalue=pvalue2, koondus=koondus, bakter_osakaal=pii1, plasmiid_osakaal1=pii2, plasmiid_osakaal2=pii3, plasmiid_osakaal3=pii4, 
  bakter_katvus=par_to_kesk(tul$par[3], keskm00), plasmiid_integr_katvus=keskm_integreeritud, plasmiid_katvus=par_to_kesk(tul2$par[3], keskm00),
  plasmiid_katvus2=par_to_kesk(tul4$par[3], keskm00), size=exp(tul$par[4]), 
  piirkond=c(1, ypiir), pr=pr, k1=k1, n1=n1, k2=k2,  n2=n2  )
class(tulem)="TestPlasmid"
return(tulem)
}

print.TestPlasmid=function(x){
  print(paste("Statistic = ", round(x$teststat,3), sep=""), quote=FALSE)
  print(paste("p-value   = ",  format(x$pvalue, digits=4), sep=""), quote=FALSE)
  print(paste("Convergence: ", x$koondus), quote=FALSE)
}

plot.TestPlasmid=function(x){
# x=ah
  katvus=0:(x$piirkond[2]*1.4)
  n=rep(0, length(katvus))
  ind=x$k1 %in% katvus
  n[x$k1[ind]+1]=x$n1[ind]
  names(n)=katvus
  plot(katvus, prop.table(n), type="p", xlim=range(katvus), col="lightblue", xlab="Coverage", ylab="distribution of k-mers", pch=20)
  ind=katvus<=x$piirkond[2] & katvus>=x$piirkond[1]
  points(katvus[ind], prop.table(n)[ind], pch=20, col="blue2")
  abline(v=x$piirkond, lty=2)

  n2=rep(0, length(katvus))
  ind=x$k2 %in% katvus
  n2[x$k2[ind]+1]=x$n2[ind]
  names(n2)=katvus
  points(katvus, prop.table(n2), xlim=range(katvus), col="pink", pch=20)
  ind=katvus<=x$piirkond[2] & katvus>=x$piirkond[1]
  points(katvus[ind], prop.table(n2)[ind], col=2, pch=20)

yA1=dnbinom(katvus,  mu=0.1,  size=x$s)
yA2=dnbinom(katvus,  mu=x$bakter_katvus*1,  size=x$s)
yA3=dnbinom(katvus,  mu=x$bakter_katvus*2,  size=x$s)
yA=yA1*x$bakter_osakaal[1]+ yA2*x$bakter_osakaal[2]+ yA3*x$bakter_osakaal[3]
lines(katvus, yA, col="darkblue", lwd=2)

yB1=dnbinom(katvus,  mu=0.1,  size=x$s)
yB2=dnbinom(katvus,  mu=x$plasmiid_katvus*1,  size=x$s)
yB3=dnbinom(katvus,  mu=x$plasmiid_katvus*2,  size=x$s)
yB=yB1*x$plasmiid_osakaal1[1]+ yB2*x$plasmiid_osakaal1[2]+ yB3*x$plasmiid_osakaal1[3]
lines(katvus, yB, col="darkred", lwd=2)

# yC1=dnbinom(katvus,  mu=0.1,  size=x$s)
# yC2=dnbinom(katvus,  mu=x$bakter_katvus*1,  size=x$s)
# yC3=dnbinom(katvus,  mu=x$bakter_katvus*2,  size=x$s)
# yC=yC1*x$plasmiid_osakaal2[1]+ yC2*x$plasmiid_osakaal2[2]+ yC3*x$plasmiid_osakaal2[3]
# lines(katvus, yC, col="red2", lwd=1)

abline(v=x$bakter_katvus, col="blue", lty=2)
abline(v=x$plasmiid_katvus, col="red", lty=2)

legend("topright", c("bacteria", "plasmid"), pch=20, col=c("blue2", "red2"), bg="white")
}

summary.TestPlasmid=function(x){
abiandmed=data.frame(source=c("Bacteria", "Plasmid, integrated", "Plasmid, unintegrated"), coverage=round(c(x$bakter_katvus, x$plasmiid_integr_katvus, x$plasmiid_katvus2),3), 
   prop_kmers_missing=round(c(x$bakter_osakaal[1], x$plasmiid_osakaal2[1], x$plasmiid_osakaal3[1]),5), prop_kmers_normal=round(c(x$bakter_osakaal[2], x$plasmiid_osakaal2[2], x$plasmiid_osakaal3[2]),5), prop_kmers_duplicated=round(c(x$bakter_osakaal[3], x$plasmiid_osakaal2[3], x$plasmiid_osakaal3[3]),5) )

print(abiandmed, row.names=FALSE)
cat("\n")
cat("Plasmid coverage multiplier:", round(x$plasmiid_katvus2/x$bakter_katvus,2))
cat("\n")
cat("\n")
cat("Coverage range used for testing:", x$piirkond[1], "...", x$piirkond[2])
cat("\n") 
cat("Allowed coverage bias:", round(x$pr*100,2), "%")
cat("\n") 
cat("\n") 
  cat("Statistic = ", round(x$teststat,3), sep="")
cat("\n") 
  cat("p-value   = ",  format(x$pvalue, digits=4), sep="")
cat("\n") 
  cat("Convergence: ", x$koondus)
  cat("\n") 
}



# Funktsioonide defineerimise l6pp *********************************************************************************************



#####################################
# Read in data and calculate values #
#####################################

bac=read.csv2(bac_file, header=T)
plasmid=read.csv2(plasmid_file, header=T)

#Launch the function

ah=testiPlasmiidi2(k1=bac[,1], n1=bac[,2], k2=plasmid[,1], n2=plasmid[,2], readi_pikkus=read_length, pr=coverage_variation, k=word)

ah$teststat
ah$pvalue
ah$koondus

if (!is.null(output_name)){
  sink(paste(output_name,".txt", sep=""))
  summary(ah)
  sink()


  png(filename=paste(output_name,".png", sep=""))
  plot(ah)
  dev.off()
}



