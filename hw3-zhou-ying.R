##Problem1
m = 100
p1 = 0.2
power_bonf = vector()
power_holm = vector()
power_hoch = vector()
power_homm = vector()
bonf = c()
holm = c()
hoch = c()
homm = c()
bonfmean = c()
holmmean = c()
hochmean = c()
hommmean = c()
for(theta in c(4,5,6,7,8,9,10)){
  for(i in c(1:100)) {
    m1 = floor(m*p1)
    sample1 = rnorm(n = m1, mean = theta, sd =1)
    sample2 = rnorm(n = m-m1, mean = 0, sd =1 )
    sample = c(sample1, sample2)#sample of some mean 0 and some mean theta
    pvalue = pnorm(sample,lower.tail = FALSE)
    pvalue = sort(pvalue)#ordered p-value
    alpha = 0.05
    bonf_alpha = alpha/m #bonferroni p-value thredhold
    bonf = c(bonf,sum(pvalue < bonf_alpha))#number of rejects by bonf method
    j = 0;
    while(pvalue[j+1]<= alpha/(m-j)){#calculate the largest i_holm
      j = j+1;}
    holm = c(holm, j)
    s = m #hochberg step up procedure
    while(pvalue[s] > alpha/(m-s+1)){s = s-1}#calculate the largest i_hoch
    hoch = c(hoch,s)
    J_star = c()
    j_star = 1
    for(j in 1:m){
      k = 1
      while((pvalue[m-j+k] > (k*alpha)/j )&& (k<=j)){#calculate J_star
        k = k+1;
      }
      if((k-1) == j){J_star = c(J_star,j)}
    }
    index2 = c()
    j_star = max(J_star)
    for(j in 1:m){
      if(pvalue[j]<= alpha/j_star){index2 = c(index2,j)}#calculate the largest i_homm
    }
    homm = c(homm, max(index2))
  }
  bonfmean = c(bonfmean, mean(bonf))
  hochmean = c(hochmean, mean(hoch))
  holmmean = c(holmmean, mean(holm))
  hommmean = c(hommmean, mean(homm))
}
bonfmean
holmmean
hochmean
hommmean
#We can see that their power from weak to strong is: bonferroni, holm, hochberg, hommel