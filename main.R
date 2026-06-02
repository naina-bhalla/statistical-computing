df<-c(1.080, 1.245, 0.111, 1.970, 2.137, 2.643, 
      1.440, 4.519, 1.914, 0.203, 3.369, 1.210, 
      0.162, 0.138, 3.060, 0.678, 1.879, 0.221, 
      1.702, 8.334, 0.878, 1.513, 1.290, 0.610, 
      2.343) 
## fort.58 dataset
## g_lambda is the function obtained by finding alpha as a function of lambda
## a closer look on log likelihood reveals that for any lambda,
## alpha -n/(sum(log(1-exp(-lambda*df)))) maximizes the log likelihood,
## so inputting this in log likelihood gives us the g_lambda function
g_lambda<-function(l)
{
  n<-length(df)
  out<-n*log(-n/(sum(log(1-exp(-l*df)))))+n*log(l)-l*sum(df)-sum(log(1-exp(-l*df)))
  return(out)
}
## first derivative of g_lambda
d_log_likelihood<-function(X,l)
{
  n<-length(X)
  A<-sum((X*exp(-l*X))/(1-exp(-l*X)))
  B<-sum(log(1-exp(-l*X)))
  out<-n/l-sum(X)+(A)*(-(n/B)-1)
  return(out)
}
## second derivative of g_lambda
double_d_log_likelihood<-function(X,l)
{
  n<-length(X)
  A<-sum((X*exp(-l*X))/(1-exp(-l*X)))
  B<-sum(log(1-exp(-l*X)))
  C<-sum(((X^2)*(exp(-l*X)))/(1-exp(-l*X))^2)
  out<-(-n/l^2)+(C)*(n/B+1)+(A^2)*(n/B^2)
  return(out)
}
##netwon raphson algorithm
NR<-function(X,l,gl,dgl,err=1e-6,max.it=1000)
{
  old_l<-l
  new_l<-l
  for(i in 1:max.it)
  {
    old_l<-new_l
    g<-gl(X,old_l)
    dg<-dgl(X,old_l)
    if(abs(dg) < 1e-12)
    {
      return(list(par=old_l,iter=i))
    }
    new_l<-old_l-g/dg
    if(new_l<0)
    {
      new_l<-old_l/2
    }
    if(abs(new_l-old_l)<err)
    {
      return(list(par=new_l,iter=i))
    }
  }
  return(list(par=old_l,iter=max.it))
}
lambda_init<-1/mean(df)
result<-NR(df,lambda_init,d_log_likelihood,double_d_log_likelihood)
lambda_hat<-result$par
iters<-result$iter
n<-length(df)
alpha_hat<-(-n/(sum(log(1-exp(-lambda_hat*df)))))

##non parametric bootstrap
set.seed(0)
B<-1000
lambda_np<-numeric(length=B)
alpha_np<-numeric(length=B)
for(i in 1:B)
{
  df_sample<-sample(df,size=n,replace=TRUE)
  res_sample<-NR(df_sample,1/mean(df_sample),d_log_likelihood,double_d_log_likelihood)
  lambda_np[i]<-res_sample$par
  alpha_np[i]<-(-n/(sum(log(1-exp(-lambda_np[i]*df_sample)))))
}
ci_95_lambda_np<-quantile(lambda_np,probs = c(0.025,0.975))
ci_95_alpha_np<-quantile(alpha_np,probs=c(0.025,0.975))

## parametric bootstrap
set.seed(0)
B<-1000
lambda_parametric<-numeric(length=B)
alpha_parametric<-numeric(length=B)

## make a transform Y=1-exp(-lambda*X), we get Y~Beta(alpha,1)
for(i in 1:B)
{
  Y<-rbeta(n,alpha_hat,1)
  df_sample_par<-(-log(1-Y)/lambda_hat)
  res_parametric<-NR(df_sample_par,1/(mean(df_sample_par)),d_log_likelihood,double_d_log_likelihood)
  lambda_parametric[i]<-res_parametric$par
  alpha_parametric[i]<-(-n/sum(log(1-exp(-lambda_parametric[i]*df_sample_par))))
}

ci_95_lambda_p<-quantile(lambda_parametric,probs = c(0.025,0.975))
ci_95_alpha_p<-quantile(alpha_parametric,probs = c(0.025,0.975))

cat("\nNewton-Raphson Optimization Results\n")
cat(sprintf("Iterations required: %d\n", iters))
cat(sprintf("Lambda MLE: %f\n", lambda_hat))
cat(sprintf("Alpha MLE:  %f\n", alpha_hat))

cat("\n95% Confidence Intervals (Non-Parametric)\n")
cat(sprintf("Lambda: [%f, %f]\n", ci_95_lambda_np[1], ci_95_lambda_np[2]))
cat(sprintf("Alpha:  [%f, %f]\n", ci_95_alpha_np[1], ci_95_alpha_np[2]))

cat("\n95% Confidence Intervals (Parametric)\n")
cat(sprintf("Lambda: [%f, %f]\n", ci_95_lambda_p[1], ci_95_lambda_p[2]))
cat(sprintf("Alpha:  [%f, %f]\n", ci_95_alpha_p[1], ci_95_alpha_p[2]))
