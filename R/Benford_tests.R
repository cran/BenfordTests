## Tests for Benford's law
## Pearson's Chi squared test statistic
chi_square_benford<-function(x=NULL,first_digits=1,pvalmethod="asymptotic",pvalsims=10000)
{
   if(!is.numeric(x)){stop("x must be numeric.")}
   
   digits<-leading_digits(x,number=first_digits)
   n<-length(digits)
   freq_of_digits<-table(c(digits,sequence_leading(first_digits)))-1
   rel_freq_of_digits<-freq_of_digits/n
   rel_freq_of_digits_H0<-pbenf(first_digits)
   
   chi_square<-n*sum((rel_freq_of_digits-rel_freq_of_digits_H0)^2/rel_freq_of_digits_H0)
   
   if(pvalmethod=="asymptotic")
   {
		pval<-1-pchisq(chi_square,df=length(sequence_leading(first_digits))-1)
   }
   if(pvalmethod=="simulate")
   {
      compute.chisq_H0<-function(n,number=1,rel_freq_of_digits_H0)
      {
         rel_freq_of_digits<-(table(c(leading_digits(rbenf(n),number=first_digits),sequence_leading(first_digits)))-1)/n;
         chi_square<-n*sum((rel_freq_of_digits-rel_freq_of_digits_H0)^2/rel_freq_of_digits_H0)
         return(chi_square)
       }
      dist_chisquareH0<-sapply(X=rep(n,pvalsims),FUN=compute.chisq_H0,number=first_digits,rel_freq_of_digits_H0)
      pval<-1-sum(dist_chisquareH0<=chi_square)/length(dist_chisquareH0)
   }
   RVAL <- list(statistic = c(Chi_square = chi_square), p.value = pval, method = "Chi-Square Test for Benford Distribution", 
                data.name = "x")
   class(RVAL) <- "htest"
   return(RVAL)
}

##Kologormov Smirnov test (EDF type)
K_S_benford<-function(x=NULL,first_digits=1,pvalmethod="simulate",pvalsims=10000)
{
   if(!is.numeric(x)){stop("x must be numeric.")}
   
   digits<-leading_digits(x,number=first_digits)
   n<-length(digits)
   freq_of_digits<-table(c(digits,sequence_leading(first_digits)))-1
   rel_freq_of_digits<-freq_of_digits/n
   rel_freq_of_digits_H0<-pbenf(first_digits)
   
   cum_sum_Ds<-cumsum(rel_freq_of_digits)-cumsum(rel_freq_of_digits_H0)
   K_S_D<-max(max(cum_sum_Ds),abs(min(cum_sum_Ds)))*sqrt(n)
   
   if(pvalmethod=="simulate")
   {
      compute.K_S_D_H0<-function(n,number=1,rel_freq_of_digits_H0)
      {
         rel_freq_of_digits<-(table(c(leading_digits(rbenf(n),number=first_digits),sequence_leading(first_digits)))-1)/n;
         cum_sum_Ds<-cumsum(rel_freq_of_digits)-cumsum(rel_freq_of_digits_H0)
         K_S_D<-max(max(cum_sum_Ds),abs(min(cum_sum_Ds)))*sqrt(n)
         return(K_S_D)
      }
      dist_K_S_D_H0<-sapply(X=rep(n,pvalsims),FUN=compute.K_S_D_H0,number=first_digits,rel_freq_of_digits_H0)
      pval<-1-sum(dist_K_S_D_H0<=K_S_D)/length(dist_K_S_D_H0)
   }
   RVAL <- list(statistic = c(K_S_D = K_S_D), p.value = pval, method = "K-S Test for Benford Distribution", 
                data.name = "x")
   class(RVAL) <- "htest"
   return(RVAL)
}

# Chebyshev Distance Test (crit values for one digit testing first by Morrow)
Chebyshev_dist_benford<-function(x=NULL,first_digits=1,pvalmethod="simulate",pvalsims=10000)
{
   if(!is.numeric(x)){stop("x must be numeric.")}
   
   digits<-leading_digits(x,number=first_digits)
   n<-length(digits)
   freq_of_digits<-table(c(digits,sequence_leading(first_digits)))-1
   rel_freq_of_digits<-freq_of_digits/n
   rel_freq_of_digits_H0<-pbenf(first_digits)
   
   m_star<-sqrt(n)*max(abs(rel_freq_of_digits-rel_freq_of_digits_H0))
   
   if(pvalmethod=="simulate")
   {
      compute.m_star_H0<-function(n,number=1,rel_freq_of_digits_H0)
      {
         rel_freq_of_digits<-(table(c(leading_digits(rbenf(n),number=first_digits),sequence_leading(first_digits)))-1)/n;
         
         m_star<-sqrt(n)*max(abs(rel_freq_of_digits-rel_freq_of_digits_H0))
         return(m_star)
      }
      dist_m_star_H0<-sapply(X=rep(n,pvalsims),FUN=compute.m_star_H0,number=first_digits,rel_freq_of_digits_H0)
      pval<-1-sum(dist_m_star_H0<=m_star)/length(dist_m_star_H0)
   }
   RVAL <- list(statistic = c(m_star = m_star), p.value = pval, method = "Tschebycheff Distance Test for Benford Distribution", 
                data.name = "x")
   class(RVAL) <- "htest"
   return(RVAL)
}

# Euclidean Distance Test(crit values for one digit testing first by Morrow)
Euclidean_dist_benford<-function(x=NULL,first_digits=1,pvalmethod="simulate",pvalsims=10000)
{
   if(!is.numeric(x)){stop("x must be numeric.")}
   
   digits<-leading_digits(x,number=first_digits)
   n<-length(digits)
   freq_of_digits<-table(c(digits,sequence_leading(first_digits)))-1
   rel_freq_of_digits<-freq_of_digits/n
   rel_freq_of_digits_H0<-pbenf(first_digits)
   
   d_star<-sqrt(n)*sqrt(sum((rel_freq_of_digits-rel_freq_of_digits_H0)^2))
   
   if(pvalmethod=="simulate")
   {
      compute.d_star_H0<-function(n,number=1,rel_freq_of_digits_H0)
      {
         rel_freq_of_digits<-(table(c(leading_digits(rbenf(n),number=first_digits),sequence_leading(first_digits)))-1)/n;
         
         d_star<-sqrt(n)*sqrt(sum((rel_freq_of_digits-rel_freq_of_digits_H0)^2))
         return(d_star)
      }
      dist_d_star_H0<-sapply(X=rep(n,pvalsims),FUN=compute.d_star_H0,number=first_digits,rel_freq_of_digits_H0)
      pval<-1-sum(dist_d_star_H0<=d_star)/length(dist_d_star_H0)
   }
   RVAL <- list(statistic = c(d_star = d_star), p.value = pval, method = "Euclidean Distance Test for Benford Distribution", 
                data.name = "x")
   class(RVAL) <- "htest"
   return(RVAL)
}

# Freedman's modification of Watsons U^2 for the Benford distribution (originally 1 digit)
Freedman_Watson_Usquare_benford<-function(x=NULL,first_digits=1,pvalmethod="simulate",pvalsims=10000)
{
   if(!is.numeric(x)){stop("x must be numeric.")}
   
   digits<-leading_digits(x,number=first_digits)
   n<-length(digits)
   freq_of_digits<-table(c(digits,sequence_leading(first_digits)))-1
   rel_freq_of_digits<-freq_of_digits/n
   rel_freq_of_digits_H0<-pbenf(first_digits)
   
   cum_sum_Ds<-cumsum(rel_freq_of_digits)-cumsum(rel_freq_of_digits_H0)
   U_square<-(n/length(rel_freq_of_digits))*(sum(cum_sum_Ds^2)-((sum(cum_sum_Ds)^2)/length(rel_freq_of_digits)))
   
   if(pvalmethod=="simulate")
   {
      compute.U_square_H0<-function(n,number=1,rel_freq_of_digits_H0)
      {
         rel_freq_of_digits<-(table(c(leading_digits(rbenf(n),number=first_digits),sequence_leading(first_digits)))-1)/n;
         
         cum_sum_Ds<-cumsum(rel_freq_of_digits)-cumsum(rel_freq_of_digits_H0)
         U_square<-(n/length(rel_freq_of_digits))*(sum(cum_sum_Ds^2)-((sum(cum_sum_Ds)^2)/length(rel_freq_of_digits)))
         return(U_square)
      }
      dist_U_square_H0<-sapply(X=rep(n,pvalsims),FUN=compute.U_square_H0,number=first_digits,rel_freq_of_digits_H0)
      pval<-1-sum(dist_U_square_H0<=U_square)/length(dist_U_square_H0)
   }
   RVAL <- list(statistic = c(U_square = U_square), p.value = pval, method = "Freedman-Watson U-squared Test for Benford Distribution", 
                data.name = "x")
   class(RVAL) <- "htest"
   return(RVAL)
}

#Normed mean deviation test for Benfords distribution first proposed as descriptive test statistic by Judge and Schechter
J_S_avg_dev_benford<-function(x=NULL,first_digits=1,pvalmethod="simulate",pvalsims=10000)
{
   if(!is.numeric(x)){stop("x must be numeric.")}
   
   digits<-leading_digits(x,number=first_digits)
   n<-length(digits)
   mu_emp<-mean(digits)
   mu_bed<-sum(sequence_leading(number=first_digits)*pbenf(first_digits))
   a_star<-abs(mu_emp-mu_bed)/(max(sequence_leading(first_digits))-mu_bed)
         
   if(pvalmethod=="simulate")
   {
      compute.a_star_H0<-function(n,number=1,mu_bed)
      {
         
         digits<-leading_digits(rbenf(n),number=number)
         n<-length(digits)
         mu_emp<-mean(digits)
         
         a_star<-abs(mu_emp-mu_bed)/(max(sequence_leading(number))-mu_bed)
         return(a_star)
      }
      dist_a_star_H0<-sapply(X=rep(n,pvalsims),FUN=compute.a_star_H0,number=first_digits,mu_bed=mu_bed)
      pval<-1-sum(dist_a_star_H0<=a_star)/length(dist_a_star_H0)
   }
   RVAL <- list(statistic = c(a_star = a_star), p.value = pval, method = "Judge-Schechter Normed Deviation Test for Benford Distribution", 
                data.name = "x")
   class(RVAL) <- "htest"
   return(RVAL)
}


# Shapiro-Francia type (correlation based) test for Benford's distribution first proposed by Joenssen (2013)
J_stat_squ_benford<-function(x=NULL,first_digits=1,method="pearson",pvalmethod="simulate",pvalsims=10000)
{
   if(!is.numeric(x)){stop("x must be numeric.")}
   
   digits<-leading_digits(x,number=first_digits)
   n<-length(digits)
   freq_of_digits<-table(c(digits,sequence_leading(first_digits)))-1
   rel_freq_of_digits<-freq_of_digits/n
   rel_freq_of_digits_H0<-pbenf(first_digits)
   
   J_stat_squ<-cor(rel_freq_of_digits,rel_freq_of_digits_H0,method=method)
   J_stat_squ<-sign(J_stat_squ)*(J_stat_squ^2)
   if(pvalmethod=="simulate")
   {
      compute.J_stat_H0<-function(n,number=1,rel_freq_of_digits_H0)
      {
         rel_freq_of_digits<-(table(c(leading_digits(rbenf(n),number=first_digits),sequence_leading(first_digits)))-1)/n;
         
         J_stat_squ<-cor(rel_freq_of_digits,rel_freq_of_digits_H0,method=method)
		 J_stat_squ<-sign(J_stat_squ)*(J_stat_squ^2)
         return(J_stat_squ)
      }
      dist_J_stat_H0<-sapply(X=rep(n,pvalsims),FUN=compute.J_stat_H0,number=first_digits,rel_freq_of_digits_H0=rel_freq_of_digits_H0)
      pval<-sum(dist_J_stat_H0<=J_stat_squ)/length(dist_J_stat_H0)
   }
   RVAL <- list(statistic = c(J_stat_squ = J_stat_squ), p.value = pval, method = "JP-square Correlation Statistic Test for Benford Distribution", 
                data.name = "x")
   class(RVAL) <- "htest"
   return(RVAL)
}

### Aditional functions provided
## returns first "number" significant digits of numerical vector x
leading_digits<-function(x=NULL, number=1)
{
   if(!is.numeric(x)){stop("x needs to be numeric.")}
   x<-abs(x)
   return(trunc((10^((floor(log10(x))*-1)+number-1))*x))   
}

##returns the sequence of all possible leading digits for "number" leading digits
#ie 1-> 1:9; 2-> 10:99; 3-> 100:999 etc.
sequence_leading<-function(number=1)
{return(seq(from=10^(number-1),to=(10^(number))-1))}

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
   benf_table<-table(sequence_leading(digits))-1
   benf_table<-benf_table+sapply(sequence_leading(number=digits),FUN=pbenf_for_seq)
   
   return(benf_table)
}

#returns a n-long sample numbers satisfying Benford's law
rbenf<-function(n)
{
   return(10^(runif(n)))
}