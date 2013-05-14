#
#    Do not delete!
#  File name		Benford_tests.R
#  Part of:		   BenfordTests (GNU R contributed package)
#  Author:			Dieter William Joenssen
#  Copyright:		Dieter William Joenssen
#  Email:			Dieter.Joenssen@TU-Ilmenau.de
#  Created:		   16 April 2013
#  Last Update: 	14 May 2013
#  Description:	R code for Package BenfordTests. Implemented functionions include following:
#                 Actual Tests:
#                 -chisq.benftest                  ~ Chi square test for Benford's law
#                 -ks.benftest                     ~ Kologormov-Smirnov test for Benford's law
#                 -mdist.benftest                  ~ Chebyshev distance based test for Benford's law
#                 -edist.benftest                  ~ Euclidean distance based test for Benford's law
#                 -usq.benftest                    ~ Freedman Watson U square test for Benford's law
#                 -meandigit.benftest              ~ mean digit test for Benford's law
#                 -jpsq.benftest                   ~ Pearson correlation test for Benford's law (removed moved ability to choose "spearman" and "kendall" for correlation)
#                 Supporting functions:
#                 -signifd                  	   ~ returns the specified number of first significant digits
#                 -signifd.seq                	   ~ sequence of all possible k first digits(i.e., k=1 -> 1:9)
#                 -qbenf                           ~ returns full cumulative probability function for Benford's distribution (first k digits)
#                 -pbenf                           ~ returns full probability function for Benford's distribution (first k digits)
#                 -rbenf                           ~ generates a random variable that satisfies Benford's law

## Tests for Benford's law
## Pearson's Chi squared test statistic
chisq.benftest<-function(x=NULL,digits=1,pvalmethod="asymptotic",pvalsims=10000)
{
#some self-explanitory error checking
   if(!is.numeric(x)){stop("x must be numeric.")}
   pvalmethod <- pmatch(pvalmethod, c("asymptotic", "simulate"))
   if (is.na(pvalmethod)){stop("invalid 'pvalmethod' argument")}
   if((length(pvalsims)!=1)){stop("'pvalsims' argument takes only single integer!")}
   if((length(digits)!=1)){stop("'digits' argument takes only single integer!")}

#reduce the data to the specified number of first digits
   first_digits<-signifd(x,digits)
  #get the amount of values that should be tested
   n<-length(first_digits)
   #the the observed frequencies of all digits
   freq_of_digits<-table(c(first_digits,signifd.seq(digits)))-1
   #calculate the relative frequencies
   rel_freq_of_digits<-freq_of_digits/n
   #get the expected frequencies under the NULL
   rel_freq_of_digits_H0<-pbenf(digits)
   #calculate the chi square test statistic
   chi_square<-n*sum((rel_freq_of_digits-rel_freq_of_digits_H0)^2/rel_freq_of_digits_H0)

   #calc pval if using the asymptotic NULL-distribution
   if(pvalmethod==1)
   {
		pval<-1-pchisq(chi_square,df=length(signifd.seq(digits))-1)
   }
   if(pvalmethod==2)#calc pval if using the simulated NULL-distribution
   {
   #wrapper function for simulating the NULL distribution
   #this is new in version 1.0.0, calculations are done in c (very fast)
      compute_chisquare_H0<-function(n=10,digits=1,pvalsims=10)
      {
	  #allocate memory
         H0_chi_square<-rep(0,pvalsims)
         H0_chi_square<- .C("compute_H0_chi_square", H0_chi_square = as.double(H0_chi_square), digits = as.integer(digits),
                            pbenf = as.double(pbenf(digits)),qbenf=as.double(qbenf(digits)),n = as.integer(n),
                            n_sim=as.integer(pvalsims))$H0_chi_square
         return(H0_chi_square)
      }
      dist_chisquareH0<-compute_chisquare_H0(n=n,digits=digits,pvalsims=pvalsims)
	  #calculate pvalue by determeninge the amount of values in the NULL-distribution that are larger than the calculated chi_square value
      pval<-1-sum(dist_chisquareH0<=chi_square)/length(dist_chisquareH0)
   }
   #make a nice S3 object of type htest
   RVAL <- list(statistic = c(chisq = chi_square), p.value = pval, method = "Chi-Square Test for Benford Distribution", 
                data.name = deparse(substitute(x)))
   class(RVAL) <- "htest"
   return(RVAL)
}

##Kologormov Smirnov test (EDF type)
ks.benftest<-function(x=NULL,digits=1,pvalmethod="simulate",pvalsims=10000)
{
#some self-explanitory error checking
   if(!is.numeric(x)){stop("x must be numeric.")}
   pvalmethod <- pmatch(pvalmethod, c("simulate"))
   if (is.na(pvalmethod)){stop("invalid 'pvalmethod' argument")}
   if((length(pvalsims)!=1)){stop("'pvalsims' argument takes only single integer!")}
   if((length(digits)!=1)){stop("'digits' argument takes only single integer!")}
   
   #reduce the data to the specified number of first digits
   first_digits<-signifd(x,digits)
     #get the amount of values that should be tested
   n<-length(first_digits)
    #the the observed (relative)frequencies of all digits
   freq_of_digits<-table(c(first_digits,signifd.seq(digits)))-1
   rel_freq_of_digits<-freq_of_digits/n
   #get the expected frequencies under the NULL
   rel_freq_of_digits_H0<-pbenf(digits)
   
   #calculate the deviations in the cumulative sums
   cum_sum_Ds<-cumsum(rel_freq_of_digits)-cumsum(rel_freq_of_digits_H0)
   #calculate the K-S-D-statistic
   K_S_D<-max(max(cum_sum_Ds),abs(min(cum_sum_Ds)))*sqrt(n)
   
   if(pvalmethod==1)#calc pval if using the simulated NULL-distribution
   {
   #wrapper function for simulating the NULL distribution
   #this is new in version 1.0.0, calculations are done in c (very fast)
      compute_KSD_H0<-function(n=10,digits=1,pvalsims=10)
      {
	  #allocate memory
         H0_KSD<-rep(0,pvalsims)
         H0_KSD<- .C("compute_H0_KSD", H0_KSD = as.double(H0_KSD), digits = as.integer(digits),
                     pbenf = as.double(pbenf(digits)),qbenf=as.double(qbenf(digits)),n = as.integer(n),
                     n_sim=as.integer(pvalsims))$H0_KSD
         return(H0_KSD)
      }
      dist_K_S_D_H0<-compute_KSD_H0(n=n,digits=digits,pvalsims=pvalsims)
	  #calculate pvalue by determeninge the amount of values in the NULL-distribution that are larger than the calculated D value
      pval<-1-sum(dist_K_S_D_H0<=K_S_D)/length(dist_K_S_D_H0)
   }
   
   #make a nice S3 object of type htest
   RVAL <- list(statistic = c(D = K_S_D), p.value = pval, method = "K-S Test for Benford Distribution", 
                data.name = deparse(substitute(x)))
   class(RVAL) <- "htest"
   return(RVAL)
}

# Chebyshev Distance Test (crit values for one digit testing first by Morrow)
mdist.benftest<-function(x=NULL,digits=1,pvalmethod="simulate",pvalsims=10000)
{
#some self-explanitory error checking
   if(!is.numeric(x)){stop("x must be numeric.")}
   pvalmethod <- pmatch(pvalmethod, c("simulate"))
   if (is.na(pvalmethod)){stop("invalid 'pvalmethod' argument")}
   if((length(pvalsims)!=1)){stop("'pvalsims' argument takes only single integer!")}
   if((length(digits)!=1)){stop("'digits' argument takes only single integer!")}
   
   #reduce the data to the specified number of first digits
   first_digits<-signifd(x,digits)
     #get the amount of values that should be tested
   n<-length(first_digits)
    #the the observed (relative)frequencies of all digits
   freq_of_digits<-table(c(first_digits,signifd.seq(digits)))-1
   rel_freq_of_digits<-freq_of_digits/n
   #get the expected frequencies under the NULL
   rel_freq_of_digits_H0<-pbenf(digits)
   
   #calculate the m_star statisitic
   #on a personal note, this isn't very different from the KSD.
   m_star<-sqrt(n)*max(abs(rel_freq_of_digits-rel_freq_of_digits_H0))
   
   if(pvalmethod==1)
   {
      #wrapper function for simulating the NULL distribution
   #this is new in version 1.0.0, calculations are done in c (very fast)
      compute_mstar_H0<-function(n=10,digits=1,pvalsims=10)
      {
	  #allocate memory
         H0_mstar<-rep(0,pvalsims)
         H0_mstar<- .C("compute_H0_mstar", H0_mstar = as.double(H0_mstar), digits = as.integer(digits),
                       pbenf = as.double(pbenf(digits)),qbenf=as.double(qbenf(digits)),n = as.integer(n),
                       n_sim=as.integer(pvalsims))$H0_mstar
         return(H0_mstar)
      }
      dist_m_star_H0<-compute_mstar_H0(n=n,digits=digits,pvalsims=pvalsims)
	  #calculate pvalue by determeninge the amount of values in the NULL-distribution that are larger than the calculated m_star value
      pval<-1-sum(dist_m_star_H0<=m_star)/length(dist_m_star_H0)
   }
   
   #make a nice S3 object of type htest
   RVAL <- list(statistic = c(m_star = m_star), p.value = pval, method = "Chebyshev Distance Test for Benford Distribution", 
                data.name = deparse(substitute(x)))
   class(RVAL) <- "htest"
   return(RVAL)
}

# Euclidean Distance Test(crit values for one digit testing first by Morrow)
edist.benftest<-function(x=NULL,digits=1,pvalmethod="simulate",pvalsims=10000)
{
#some self-explanitory error checking
   if(!is.numeric(x)){stop("x must be numeric.")}
   pvalmethod <- pmatch(pvalmethod, c("simulate"))
   if (is.na(pvalmethod)){stop("invalid 'pvalmethod' argument")}
   if((length(pvalsims)!=1)){stop("'pvalsims' argument takes only single integer!")}
   if((length(digits)!=1)){stop("'digits' argument takes only single integer!")}
   
   #reduce the data to the specified number of first digits
   first_digits<-signifd(x,digits)
     #get the amount of values that should be tested
   n<-length(first_digits)
    #the the observed (relative)frequencies of all digits
   freq_of_digits<-table(c(first_digits,signifd.seq(digits)))-1
   rel_freq_of_digits<-freq_of_digits/n
   #get the expected frequencies under the NULL
   rel_freq_of_digits_H0<-pbenf(digits)
   
   #calculate the m_star statisitic
   #on a personal note, this isn't very different from the m_star statistic.
   d_star<-sqrt(n)*sqrt(sum((rel_freq_of_digits-rel_freq_of_digits_H0)^2))
   
   if(pvalmethod==1)
   {
   #wrapper function for simulating the NULL distribution
   #this is new in version 1.0.0, calculations are done in c (very fast)
      compute_dstar_H0<-function(n=10,digits=1,pvalsims=10)
      {
	  #allocate memory
         H0_dstar<-rep(0,pvalsims)
         H0_dstar<- .C("compute_H0_dstar", H0_dstar = as.double(H0_dstar), digits = as.integer(digits),
                       pbenf = as.double(pbenf(digits)),qbenf=as.double(qbenf(digits)),n = as.integer(n),
                       n_sim=as.integer(pvalsims))$H0_dstar
         return(H0_dstar)
      }
      dist_d_star_H0<-compute_dstar_H0(n=n,digits=digits,pvalsims=pvalsims)
	  #calculate pvalue by determeninge the amount of values in the NULL-distribution that are larger than the calculated d_star value
      pval<-1-sum(dist_d_star_H0<=d_star)/length(dist_d_star_H0)
   }
   
   #make a nice S3 object of type htest
   RVAL <- list(statistic = c(d_star = d_star), p.value = pval, method = "Euclidean Distance Test for Benford Distribution", 
                data.name = deparse(substitute(x)))
   class(RVAL) <- "htest"
   return(RVAL)
}

# Freedman's modification of Watsons U^2 for the Benford distribution (originally 1 digit)
usq.benftest<-function(x=NULL,digits=1,pvalmethod="simulate",pvalsims=10000)
{
 #some self-explanitory error checking
   if(!is.numeric(x)){stop("x must be numeric.")}
   pvalmethod <- pmatch(pvalmethod, c("simulate"))
   if (is.na(pvalmethod)){stop("invalid 'pvalmethod' argument")}
   if((length(pvalsims)!=1)){stop("'pvalsims' argument takes only single integer!")}
   if((length(digits)!=1)){stop("'digits' argument takes only single integer!")}
   
   #reduce the data to the specified number of first digits
   first_digits<-signifd(x,digits)
     #get the amount of values that should be tested
   n<-length(first_digits)
    #the the observed (relative)frequencies of all digits
   freq_of_digits<-table(c(first_digits,signifd.seq(digits)))-1
   rel_freq_of_digits<-freq_of_digits/n
   #get the expected frequencies under the NULL
   rel_freq_of_digits_H0<-pbenf(digits)
   
   #calculate deviations betwen the cumulative sums
   cum_sum_Ds<-cumsum(rel_freq_of_digits)-cumsum(rel_freq_of_digits_H0)
   #calculate the U^2 test statistic
   U_square<-(n/length(rel_freq_of_digits))*(sum(cum_sum_Ds^2)-((sum(cum_sum_Ds)^2)/length(rel_freq_of_digits)))
   
   if(pvalmethod==1)
   {
      #wrapper function for simulating the NULL distribution
   #this is new in version 1.0.0, calculations are done in c (very fast)
      compute_U_square_H0<-function(n=10,digits=1,pvalsims=10)
      {
	  #allocate memory
         H0_U_square<-rep(0,pvalsims)
         H0_U_square<- .C("compute_H0_U_square", H0_U_square = as.double(H0_U_square), digits = as.integer(digits),
                          pbenf = as.double(pbenf(digits)),qbenf=as.double(qbenf(digits)),n = as.integer(n),
                          n_sim=as.integer(pvalsims))$H0_U_square
         return(H0_U_square)
      }
      dist_U_square_H0<-compute_U_square_H0(n=n,digits=digits,pvalsims=pvalsims)
	  #calculate pvalue by determeninge the amount of values in the NULL-distribution that are larger than the calculated U_square value
      pval<-1-sum(dist_U_square_H0<=U_square)/length(dist_U_square_H0)
   }
   
   #make a nice S3 object of type htest
   RVAL <- list(statistic = c(U_square = U_square), p.value = pval, method = "Freedman-Watson U-squared Test for Benford Distribution", 
                data.name = deparse(substitute(x)))
   class(RVAL) <- "htest"
   return(RVAL)
}

#Normed mean deviation test for Benfords distribution first proposed as descriptive test statistic by Judge and Schechter
meandigit.benftest<-function(x=NULL,digits=1,pvalmethod="simulate",pvalsims=10000)
{
#some self-explanitory error checking
   if(!is.numeric(x)){stop("x must be numeric.")}
   pvalmethod <- pmatch(pvalmethod, c("simulate"))
   if (is.na(pvalmethod)){stop("invalid 'pvalmethod' argument")}
   if((length(pvalsims)!=1)){stop("'pvalsims' argument takes only single integer!")}
   if((length(digits)!=1)){stop("'digits' argument takes only single integer!")}
   
   #get specified number of first digits and number of numbers
   first_digits<-signifd(x,digits)
   n<-length(first_digits)
   #get empirical mean digit value
   mu_emp<-mean(first_digits)
   #get expected mean digit value under NULL
   mu_bed<-sum(signifd.seq(digits)*pbenf(digits))
   #normalize to get a_star
   a_star<-abs(mu_emp-mu_bed)/(max(signifd.seq(digits))-mu_bed)
         
   if(pvalmethod==1)
   {
         #wrapper function for simulating the NULL distribution
   #this is new in version 1.0.0, calculations are done in c (very fast)
      compute_astar_H0<-function(n=10,digits=1,pvalsims=10)
      {
	  #allocate memory
         H0_astar<-rep(0,pvalsims)
         H0_astar<- .C("compute_H0_astar", H0_astar = as.double(H0_astar), digits = as.integer(digits),
                       pbenf = as.double(pbenf(digits)),qbenf=as.double(qbenf(digits)),n = as.integer(n),
                       n_sim=as.integer(pvalsims))$H0_astar
         return(H0_astar)
      }
      dist_a_star_H0<-compute_astar_H0(n=n,digits=digits,pvalsims=pvalsims)
	  #calculate pvalue by determeninge the amount of values in the NULL-distribution that are larger than the calculated a_star value
      pval<-1-sum(dist_a_star_H0<=a_star)/length(dist_a_star_H0)
	  #this is a two-sided test p_value must be adjusted properly
	  if(pval>.5){pval<- (1- pval)*2}
	  else{pval<- pval*2}
   }
   
   #make a nice S3 object of type htest
   RVAL <- list(statistic = c(a_star = a_star), p.value = pval, method = "Judge-Schechter Normed Deviation Test for Benford Distribution", 
                data.name = deparse(substitute(x)))
   class(RVAL) <- "htest"
   return(RVAL)
}


# Shapiro-Francia type (correlation based) test for Benford's distribution first proposed by Joenssen (2013)
jpsq.benftest<-function(x=NULL,digits=1,pvalmethod="simulate",pvalsims=10000)
{
 #some self-explanitory error checking
   if(!is.numeric(x)){stop("x must be numeric.")}
   pvalmethod <- pmatch(pvalmethod, c("simulate"))
   if (is.na(pvalmethod)){stop("invalid 'pvalmethod' argument")}
   if((length(pvalsims)!=1)){stop("'pvalsims' argument takes only single integer!")}
   if((length(digits)!=1)){stop("'digits' argument takes only single integer!")}
   
   #reduce the data to the specified number of first digits
   first_digits<-signifd(x,digits)
     #get the amount of values that should be tested
   n<-length(first_digits)
    #the the observed (relative)frequencies of all digits
   freq_of_digits<-table(c(first_digits,signifd.seq(digits)))-1
   rel_freq_of_digits<-freq_of_digits/n
   #get the expected frequencies under the NULL
   rel_freq_of_digits_H0<-pbenf(digits)
   
   #calculate the Jstat statistic
   J_stat_squ<-cor(rel_freq_of_digits,rel_freq_of_digits_H0)
   #square the J_stat_statistic and adjust the sign
   J_stat_squ<-sign(J_stat_squ)*(J_stat_squ^2)
   if(pvalmethod==1)
   {
   #wrapper function for simulating the NULL distribution
   #this is new in version 1.0.0, calculations are done in c (very fast)
      compute_J_stat_H0<-function(n=10,digits=1,pvalsims=10)
      {
	  #allocate memory
         H0_J_stat<-rep(0,pvalsims)
         H0_J_stat<- .C("compute_H0_J_stat", H0_J_stat = as.double(H0_J_stat), digits = as.integer(digits),
                        pbenf = as.double(pbenf(digits)),qbenf=as.double(qbenf(digits)),n = as.integer(n),
                        n_sim=as.integer(pvalsims))$H0_J_stat
         return(H0_J_stat)
      }
      dist_J_stat_H0<- compute_J_stat_H0(n=n,digits=digits,pvalsims=pvalsims)
	  #calculate pvalue by determeninge the amount of values in the NULL-distribution that are larger than the calculated J_stat_square value
      pval<-sum(dist_J_stat_H0<=J_stat_squ)/length(dist_J_stat_H0)
   }
   
   #make a nice S3 object of type htest
   RVAL <- list(statistic = c(J_stat_squ = J_stat_squ), p.value = pval, method = "JP-Square Correlation Statistic Test for Benford Distribution", 
                data.name = deparse(substitute(x)))
   class(RVAL) <- "htest"
   return(RVAL)
}

### Aditional functions provided
## returns first "digits" significant digits of numerical vector x
signifd<-function(x=NULL, digits=1)
{
#some self-explanitory error checking
   if(!is.numeric(x)){stop("x needs to be numeric.")}
   #calculate the first significant digits
   x<-abs(x)
   return(trunc((10^((floor(log10(x))*-1)+digits-1))*x))   
}

##returns the sequence of all possible leading digits for "digits" leading digits
#ie 1-> 1:9; 2-> 10:99; 3-> 100:999 etc.
signifd.seq<-function(digits=1)
{return(seq(from=10^(digits-1),to=(10^(digits))-1))}

# returns complete cumulative distribution function of Benford distribution for the given amount of significant digits
qbenf<-function(digits=1)
{
	return(cumsum(pbenf(digits)))
}

# returns complete probability distribution function of Benford distribution for the given amount of significant digits
pbenf<-function(digits=1)
{
   pbenf_for_seq<-function(leaddigit=10)
   {
      return(log10(1+(1/leaddigit)))
   }
   benf_table<-table(signifd.seq(digits))-1
   benf_table<-benf_table+sapply(signifd.seq(digits),FUN=pbenf_for_seq)
   
   return(benf_table)
}

#returns a n-long sample numbers satisfying Benford's law
rbenf<-function(n)
{
   return(10^(runif(n)))
}



