#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>
#include <RcppParallel.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]


using namespace Rcpp ;
using namespace RcppParallel;


//########################################################################################################################################
// function to generate normal random variable
// [[Rcpp::export]]
double rnorm(double a, double b) { // a: mean, b: s.d.
	double c=a+b*sum(rnorm(1));
    return c;
}

//########################################################################################################################################
// function to compute the nb distribution
// [[Rcpp::export]]
double dnb(int n, int m, double c) { 
	
double d=R::dpois(n,(m+1)/c,1);
return d;
}

//########################################################################################################################################
// function to for inverse of logit
// [[Rcpp::export]]
double logit(double a) { 
    return log(a/(1-a));
}

//########################################################################################################################################
// function to for inverse of logit
// [[Rcpp::export]]
double invlogit(double a) { 
    return exp(a)/(1+exp(a));
}


//########################################################################################################################################
// function to generate binomial random number
// [[Rcpp::export]]
int gen_binom(double p){
double cut=(double)rand()/(RAND_MAX);
int out=0;
if (cut<p){
out=1;
}
return out;
}

//########################################################################################################################################
// function to general multiviariable normal given sigma
// [[Rcpp::export]]
NumericVector rmnorm(arma::mat sigma) {
int ncols=sigma.n_cols;
arma::rowvec c=arma::randn(1,ncols);
arma::rowvec a=c*arma::chol(sigma);   
NumericVector b=NumericVector(a.begin(),a.end());   
return b;
}




//########################################################################################################################################
// function for simulation
// [[Rcpp::export]]
List sim_data(IntegerMatrix data, 
IntegerMatrix contactmatrix,	
NumericMatrix data2,
NumericVector para, 
NumericVector para2,
NumericVector SI_vector,	
NumericVector inc_vector,
NumericVector missing){ 

int b0;
int b1;
int b2;
int b3;
int b4;
double cut=0;
// clone the input data for simulation
IntegerMatrix data1(clone(data));
NumericMatrix record(data1.nrow(),8);

// generate the age and sex
for (b1=data1.nrow()-1;b1>=0;--b1){	
// medppl
data1(b1,19)=gen_binom(para2[31]);
// missing sex
//if (data1(b1,8)==-1){
data1(b1,8)=gen_binom(para2[2*(data1(b1,7)-1)+1]);
//}
// missing age
//if (data1(b1,9)==-1){
cut=(double)rand()/(RAND_MAX);
for (b2=7;b2>=0;--b2){
record(b1,b2)=para2[6+8*(data1(b1,7)-1)+b2];
cut-=para2[6+8*(data1(b1,7)-1)+b2];
if (cut<0){
data1(b1,9)=b2;
break;
}
//}
}
// remove those information 
// non-primary case
if (data1(b1,1)==0){
// update exposure bound
data1(b1,12)=data1(b1,15);
data1(b1,13)=data1(b1,16);
data1(b1,10)=0;
data1(b1,11)=-1;
data1(b1,14)=-1;
}
// primary case
else{
cut=(double)rand()/(RAND_MAX);
if (data1(b1,11)==-1){
data1(b1,14)=data1(b1,12)+rand()%(data1(b1,13)-data1(b1,12)+1);
// inc
for (b2=13;b2>=0;--b2){
cut-=inc_vector[b2];
if (cut<0){
data1(b1,11)=data1(b1,14)+b2+1;	
break;
}
}
}
else{
for (b2=13;b2>=0;--b2){
cut-=inc_vector[b2];
if (cut<0){
data1(b1,14)=data1(b1,11)-b2-1;	
break;
}
}
}
}
}

// hazard matrix (it is escape probability matrix)
NumericMatrix hazard(data1.nrow(),150); 	
// create sus vector
NumericVector sus(data1.nrow());
NumericMatrix effrisk(data1.nrow(),2);
for (b1=data1.nrow()-1;b1>=0;--b1){
if (data1(b1,1)==0){
// define sus	
// sex
sus(b1)=para[8]*(data1(b1,8)==1);
// age
sus(b1)+=para[9]*(data1(b1,9)/2==0);
sus(b1)+=para[10]*(data1(b1,9)/2==2);
sus(b1)+=para[11]*(data1(b1,9)/2==3);

sus(b1)+=para[12]*(data1(b1,7)==2);
sus(b1)+=para[13]*(data1(b1,7)==3);
sus(b1)+=para[14]*(data1(b1,18)==1);
sus(b1)+=para[15]*(data1(b1,18)==3);
sus(b1)+=para[16]*(data1(b1,19)==1);
// contactype
effrisk(b1,1)=invlogit(logit(para[data1(b1,3)])+sus(b1));
effrisk(b1,0)=invlogit(logit(para[0])+sus(b1));
// add a small hazard in the expousre period
for (b3=140;b3>=1;--b3){
if ((b3>=data1(b1,15))&&(b3<=data1(b1,16))){
hazard(b1,b3-1)=(1-effrisk(b1,0));
}	
}	
// add the force for primary
// screen all other people and add the force the the vector
for (b4=contactmatrix.ncol()-1;b4>=0;--b4){
if (contactmatrix(b1,b4)!=-1){	
b2=contactmatrix(b1,b4)-1;
if (data1(b2,10)==1){	
// origin is symptom onset
for (b3=data1(b2,11)-5;b3<=data1(b2,11)-5+SI_vector.length()-1;++b3){
if ((b3>=data1(b1,15))&&(b3<=data1(b1,16))&&(b3>data1(b2,14))){
double prerisk=invlogit(data2(b2,0)+logit(para[data1(b1,3)])+sus[b1]+para[7]*(b3<data1(b2,11)));
hazard(b1,b3-1)*=1-SI_vector[b3-data1(b2,11)+5]*prerisk;
}
}
}
}
}
}
}
// b2 is time here
for (b0=1;b0<=140;++b0){
for (b1=data1.nrow()-1;b1>=0;--b1){
if ((data1(b1,1)==0)&&(data1(b1,10)==0)){	
double temp_hazard=hazard(b1,b0-1); //*exp(sus(b1));
if (temp_hazard>0){
if (gen_binom(1-temp_hazard)){	
// inf date
data1(b1,10)=1;
data1(b1,14)=b0;
// onset date
cut=(double)rand()/(RAND_MAX);
for (b2=13;b2>=0;--b2){
cut-=inc_vector[b2];
if (cut<0){
data1(b1,11)=data1(b1,14)+b2+1;	
// update upperbound
if (data1(b1,11)-1<data1(b1,13)){
data1(b1,13)=data1(b1,11)-1;
}
break;
}
}
}
}
}
}
// update the hazard profile
for (b1=data1.nrow()-1;b1>=0;--b1){
for (b4=contactmatrix.ncol()-1;b4>=0;--b4){
if (contactmatrix(b1,b4)!=-1){	
b2=contactmatrix(b1,b4)-1;
if ((data1(b2,10)==1)&&(data1(b2,14)==b0)&&(data1(b2,1)==0)){	
// origin is symptom onset
for (b3=data1(b2,11)-5;b3<=data1(b2,11)-5+SI_vector.length()-1;++b3){
if ((b3>=data1(b1,15))&&(b3<=data1(b1,16))&&(b3>data1(b2,14))){
double prerisk=invlogit(data2(b2,0)+logit(para[data1(b1,3)])+sus[b1]+para[7]*(b3<data1(b2,11)));
hazard(b1,b3-1)*=1-SI_vector[b3-data1(b2,11)+5]*prerisk;	
}
}
}
}
}
}

}

// add missing
for (b1=data1.nrow()-1;b1>=0;--b1){
if (missing(0)){
if (data(b1,8)==-1){
data1(b1,8)=-1;
}
if (data(b1,9)==-1){
data1(b1,9)=-1;
}	
}
if (missing(1)){
if (data1(b1,10)==1){
if (data1(b1,12)!=data1(b1,13)){
data1(b1,14)=-1;	
}
if (gen_binom(0.02)){
data1(b1,11)=-1;
}
// inf1
if (gen_binom(0.16)){
data1(b1,18)=-1;
}
// inf2
if (gen_binom(0.08)){
data1(b1,19)=-1;
}

}	
}
}

return List::create(_[""]=data1,
	_[""]=hazard,
	_[""]=sus,
	_[""]=effrisk,
	_[""]=record);
} 



//########################################################################################################################################
//########################################################################################################################################
// the main body of the parallel for Digraphlikelihood
struct LogLik:public Worker{
// source vector
RMatrix<double> out;
RMatrix<int> data1;
RMatrix<int> contactmatrix;
RMatrix<double> data2;
RVector<double> para;
RVector<double> para2;
RVector<double> SI_vector;
RVector<double> inc_vector;
RMatrix<double> hazard;
RVector<double> sus;
RMatrix<double> effrisk;
// destination vector
// initialize with source and destination
LogLik(NumericMatrix out,
IntegerMatrix data1,
IntegerMatrix contactmatrix,
NumericMatrix data2,
NumericVector para,
NumericVector para2,
NumericVector SI_vector,
NumericVector inc_vector,
NumericMatrix hazard,
NumericVector sus,
NumericMatrix effrisk) 
:out(out),data1(data1),contactmatrix(contactmatrix),data2(data2),para(para),para2(para2),SI_vector(SI_vector),inc_vector(inc_vector),hazard(hazard),sus(sus),effrisk(effrisk){}
void operator()(std::size_t begin, std::size_t end) {


// section to write the parallel version
// functor (pass input and output matrixes)
for (unsigned int b1=begin;b1<end;++b1){
	
if (data1(b1,1)==0){

int b2;
int b3;
int b4;

// first compute his susceptibility
// sex
sus[b1]=para[8]*(data1(b1,8)==1);
// age
sus[b1]+=para[9]*(data1(b1,9)/2==0);	
sus[b1]+=para[10]*(data1(b1,9)/2==2);	
sus[b1]+=para[11]*(data1(b1,9)/2==3);	

sus[b1]+=para[12]*(data1(b1,7)==2);
sus[b1]+=para[13]*(data1(b1,7)==3);
sus[b1]+=para[14]*(data1(b1,18)==1);
sus[b1]+=para[15]*(data1(b1,18)==3);
sus[b1]+=para[16]*(data1(b1,19)==1);

// contactype
effrisk(b1,1)=invlogit(logit(para[data1(b1,3)])+sus[b1]);
effrisk(b1,0)=invlogit(logit(para[0])+sus[b1]);

// add a small hazard in the expousre period

for (b3=149;b3>=1;--b3){
if ((b3>=data1(b1,12))&&(b3<=data1(b1,17))){
hazard(b1,b3-1)=1-effrisk(b1,0);
}
}

// screen all other people and add the force the the vector
for (b4=contactmatrix.ncol()-1;b4>=0;--b4){
if (contactmatrix(b1,b4)!=-1){	
b2=contactmatrix(b1,b4)-1;
if (data1(b2,10)==1){	
// origin is symptom onset
for (b3=data1(b2,11)-5;b3<=data1(b2,11)-5+SI_vector.length()-1;++b3){
if ((b3>=data1(b1,12))&&(b3<=data1(b1,13))&&(b3>data1(b2,14))){
double prerisk=invlogit(data2(b2,0)+logit(para[data1(b1,3)])+sus[b1]+para[7]*(b3<data1(b2,11)));
hazard(b1,b3-1)*=1-SI_vector[b3-data1(b2,11)+5]*prerisk;
}
}
}
}
}

// output the likelihood
if (data1(b1,10)==0){
for (b2=data1(b1,12);b2<=data1(b1,17);++b2){
out(b1,b2-1)=log(hazard(b1,b2-1)); //*exp(sus);
}
}

if (data1(b1,10)==1){
for (b2=data1(b1,12);b2<data1(b1,14);++b2){
out(b1,b2-1)=log(hazard(b1,b2-1)); //*exp(sus);
}
out(b1,data1(b1,14)-1)=log(1-hazard(b1,data1(b1,14)-1));
}
}

// incubation period likelihood
if (data1(b1,10)==1){
// add the incubation period at column 1
int inclength=data1(b1,11)-data1(b1,14);
if ((inclength<inc_vector.length())&&(inclength>0)){
out(b1,0)+=log(inc_vector[inclength-1]);	
}
else{
out(b1,0)+=log(0.000000000001);
}
}

// age sex likelihood

out(b1,1)+=log(para2[2*(data1(b1,7)-1)+data1(b1,8)]);
out(b1,2)+=log(para2[6+8*(data1(b1,7)-1)+data1(b1,9)]);
out(b1,3)+=log(para2[30+data1(b1,19)]);

// random effect likelihood
if (contactmatrix(b1,0)==-1){
out(b1,4)+=R::dnorm(data2(b1,0),0.0,para[17],1);	
}

}
}
};

//########################################################################################################################################
//function to compute likelihood for infection and symptom
// [[Rcpp::export]]
NumericMatrix loglik(IntegerMatrix data1,
IntegerMatrix contactmatrix,
NumericMatrix data2,	
NumericVector para,
NumericVector para2,
NumericVector SI_vector,
NumericVector inc_vector){
// data2 is random effect, should only non-zero for primary cases
// record the likelihood value
NumericMatrix hazard(data1.nrow(),150); 	
NumericVector sus(data1.nrow());
NumericMatrix effrisk(data1.nrow(),2);
NumericMatrix out(data1.nrow(),150); 

// create a vector to compute
LogLik loglik(out,data1,contactmatrix,data2,para,para2,SI_vector,inc_vector,hazard,sus,effrisk);
// call parallelFor to do the work
parallelFor(0,data1.nrow(),loglik);

return out;
}

//########################################################################################################################################
//function to compute likelihood for infection and symptom
// [[Rcpp::export]]
List loglik2(IntegerMatrix data1,
IntegerMatrix contactmatrix,
NumericMatrix data2,	
NumericVector para,
NumericVector para2,
NumericVector SI_vector,
NumericVector inc_vector){

// record the likelihood value
NumericMatrix hazard(data1.nrow(),150); 	
NumericVector sus(data1.nrow());
NumericMatrix effrisk(data1.nrow(),2);
NumericMatrix out(data1.nrow(),150); 

// create a vector to compute
LogLik loglik(out,data1,contactmatrix,data2,para,para2,SI_vector,inc_vector,hazard,sus,effrisk);
// call parallelFor to do the work
parallelFor(0,data1.nrow(),loglik);

return List::create(_[""]=out,
	_[""]=hazard,
	_[""]=sus,
	_[""]=effrisk);
}


//########################################################################################################################################
//function to compute the prior likelihood 
// [[Rcpp::export]]
double prior_loglik(NumericVector para){
// check if the para are within their possible range
NumericVector out(para.length());
int b1;
for (b1=5;b1>=0;--b1){
out(b1)=R::dunif(para(b1),0.000000000000000001,0.99,1);
}
for (b1=16;b1>=6;--b1){
out(b1)=R::dunif(para(b1),-10.0,10.0,1);
}
out(17)=R::dgamma(para(17),0.001,1/0.001,1);  
double output=sum(out);
// if the prior is outside the parameter space
if (output< -9999999){
output=-9999999;
}
return output;
}



//##############################################################################################################################################
//##############################################################################################################################################
// function for mcmc
// [[Rcpp::export]]
List mcmc(IntegerMatrix data1, 
IntegerMatrix contactmatrix,	
NumericVector int_para, 
NumericVector int_para2,
NumericVector SI_vector,	
NumericVector inc_vector,
int mcmc_n,             
NumericVector move,    
NumericVector sigma){            

// create the vector for use
int b0;
int b1;
int b2;
int b3;
int b4;
int moveindex;

// need to set number of parameter here
NumericMatrix p_para(mcmc_n,int_para.length());
NumericMatrix p_para_r(mcmc_n,sum(move));
p_para(0,_)=int_para;
moveindex=sum(move)-1;
for (b1=int_para.length()-1;b1>=0;--b1){
if (move(b1)){
p_para_r(0,moveindex)=p_para(0,b1);
--moveindex;
}	
}

// create random effect   
NumericMatrix data2(data1.nrow(),1);
for (b1=data2.nrow()-1;b1>=0;--b1){
if (contactmatrix(b1,0)==-1){
data2(b1,0)=R::rnorm(0.0,int_para[17]);	
}	
}

// for the second parameter vectors
NumericMatrix baselinecount(1,int_para2.length());
// count the real count here because those would not change
for (b1=data1.nrow()-1;b1>=0;--b1){
if (data1(b1,1)!=1){	
// male
if (data1(b1,8)!=-1){
++baselinecount(0,2*(data1(b1,7)-1)+data1(b1,8));
}
// age group
if (data1(b1,9)!=-1){
++baselinecount(0,6+8*(data1(b1,7)-1)+data1(b1,9));
}
// medppl
if ((data1(b1,19)!=-1)){
++baselinecount(0,30+data1(b1,19));
}
}
}


NumericMatrix p_para2(mcmc_n,int_para2.length());
p_para2(0,_)=int_para2;


// create augmented data
IntegerMatrix data11(clone(data1));

for (b1=data11.nrow()-1;b1>=0;--b1){
// impute the infection time	
if ((data1(b1,14)==-1)&&(data1(b1,10)==1)){	
data11(b1,14)=data11(b1,13);	
}
// imput sympotm onset
if ((data1(b1,11)==-1)&&(data1(b1,10)==1)){	
data11(b1,11)=data11(b1,14)+1;	
}
// impute age
if (data1(b1,9)==-1){
data11(b1,9)=0;	
}
// impute sex
if (data1(b1,8)==-1){
data11(b1,8)=0;
}
// impute inf1 and inf2
if (data1(b1,19)==-1){
data11(b1,19)=1;	
}
}


// matrix to record LL
NumericMatrix LL1(mcmc_n,3);

//####################################################################################################################################
// compute likelihood

NumericMatrix loglik1=loglik(data11,contactmatrix,data2,p_para(0,_),p_para2(0,_),SI_vector,inc_vector);
NumericMatrix loglik1pro;

LL1(0,0)=sum(loglik1);

NumericVector temploglik(1);
NumericVector newloglik(1);
temploglik(0)=LL1(0,0)+prior_loglik(p_para(0,_));

double loglikeratio;
double accept_pro;
NumericVector pro_para(int_para.length());

// record rm
NumericMatrix rmrecord(mcmc_n,data1.nrow());
for (b1=data1.nrow()-1;b1>=0;--b1){
if (contactmatrix(b1,0)==-1){
rmrecord(0,b1)=data2(b1,0);	
}	
}

//####################################################################################################################################
// main mcmc step
NumericVector acceptrate(int_para.length());
NumericMatrix imputeacceptrate(data11.nrow(),6);
IntegerMatrix infection_record(mcmc_n,sum(data1(_,10)));

//####################################################################################################################################
for (b0=1;b0<mcmc_n;++b0){


// after 500 step, then set the sigma to be the empirical sigma
if ((b0>500)&&(b0%100==0)){
for (b1=int_para.length()-1;b1>=0;--b1){
if (move(b1)){	
NumericVector temp1(b0-1);
for (b2=b0-2;b2>=0;--b2){
temp1(b2)=p_para(b2,b1);	
}
sigma(b1)=sd(temp1);
// tuning
if (acceptrate(b1)<0.1){
sigma(b1)*=0.5;
}	
if ((acceptrate(b1)<0.15)&(acceptrate(b1)>0.1)){
sigma(b1)*=0.8;
}
if ((acceptrate(b1)<0.2)&(acceptrate(b1)>0.15)){
sigma(b1)*=0.95;
}
if ((acceptrate(b1)<0.4)&(acceptrate(b1)>0.3)){
sigma(b1)*=1.05;
}
if ((acceptrate(b1)<0.9)&(acceptrate(b1)>0.4)){
sigma(b1)*=1.2;
}
if (acceptrate(b1)>0.9){
sigma(b1)*=2;
}
}
}

}

// metorpolis-hasing update on parameter

for (b1=0;b1<int_para.length();++b1){
if (move(b1)){
pro_para=p_para(b0-1,_);
for (b2=b1-1;b2>=0;--b2){
pro_para(b2)=p_para(b0,b2);	
}
pro_para(b1)+=rnorm(0.0,sigma(b1));
newloglik(0)=prior_loglik(pro_para);
if (newloglik(0)> -9999999){
loglik1pro=loglik(data11,contactmatrix,data2,pro_para,p_para2(b0-1,_),SI_vector,inc_vector);
newloglik(0)+=sum(loglik1pro);
loglikeratio=newloglik(0)-temploglik(0);
accept_pro=pow(exp(1),loglikeratio);
}
else{
accept_pro=0;	
}
if(gen_binom(accept_pro)){
loglik1=clone(loglik1pro);		
temploglik(0)=newloglik(0);
p_para(b0,b1)=pro_para(b1);
acceptrate(b1)*=(b0-1);
acceptrate(b1)+=1;
acceptrate(b1)/=b0;
}
else{
p_para(b0,b1)=p_para(b0-1,b1);
acceptrate(b1)*=(b0-1);
acceptrate(b1)/=b0;
}
}
else {
p_para(b0,b1)=p_para(b0-1,b1);
}
}

LL1(b0,0)=temploglik(0)-prior_loglik(p_para(b0,_));

// gibbs sampling for age and sex distribution
// here to update para2
NumericVector tempcount(7);
for (b1=baselinecount.ncol()-1;b1>=0;--b1){
int para3index=1;
if (b1<=1){
para3index=0;	
}
if ((b1>=4)&&(b1<=5)){
para3index=2;	
}
if ((b1>=6)&&(b1<=13)){
para3index=3;
}
if ((b1>=14)&&(b1<=21)){
para3index=4;	
}
if ((b1>=22)&&(b1<=29)){
para3index=5;
}
if ((b1>=30)&&(b1<=31)){
para3index=6;
}
p_para2(b0,b1)=R::rgamma(baselinecount(0,b1)+1.0,1.0);
tempcount(para3index)+=p_para2(b0,b1);
}
for (b1=baselinecount.ncol()-1;b1>=0;--b1){
int para3index=1;
if (b1<=1){
para3index=0;	
}
if ((b1>=4)&&(b1<=5)){
para3index=2;	
}
if ((b1>=6)&&(b1<=13)){
para3index=3;
}
if ((b1>=14)&&(b1<=21)){
para3index=4;	
}
if ((b1>=22)&&(b1<=29)){
para3index=5;
}	
if ((b1>=30)&&(b1<=31)){
para3index=6;
}
p_para2(b0,b1)/=tempcount(para3index);
}

// update likelihood after update para2

loglik1=loglik(data11,contactmatrix,data2,p_para(b0,_),p_para2(b0,_),SI_vector,inc_vector);

LL1(b0,0)=sum(loglik1);
temploglik(0)=LL1(b0,0)+prior_loglik(p_para(b0,_));

// sample missing infomration
// sample age and sex
for (b1=data1.nrow()-1;b1>=0;--b1){
if (data1(b1,1)!=1){
int update=0;
if ((data1(b1,8)==-1)||(data1(b1,9)==-1)||(data1(b1,19)==-1)){
update=gen_binom(0.1);	
}
if (update){
accept_pro=1;
IntegerMatrix data11_pro(clone(data11));	
// sample sex
if (data1(b1,8)==-1){
data11_pro(b1,8)=1-data11(b1,8);
}
// sample age
if (data1(b1,9)==-1){
data11_pro(b1,9)=rand()%8;
}
// sample medppl
if (data1(b1,19)==-1){
data11_pro(b1,19)=1-data11(b1,19);
}
newloglik(0)=prior_loglik(p_para(b0,_));
if (newloglik(0)> -9999999){
loglik1pro=loglik(data11_pro,contactmatrix,data2,p_para(b0,_),p_para2(b0,_),SI_vector,inc_vector);
newloglik(0)+=sum(loglik1pro);
loglikeratio=newloglik(0)-temploglik(0);
accept_pro=pow(exp(1),loglikeratio);
}
if(gen_binom(accept_pro)){
loglik1=clone(loglik1pro);		
temploglik(0)=newloglik(0);
data11=clone(data11_pro);
imputeacceptrate(b1,0)*=(b0-1);
imputeacceptrate(b1,0)+=1;
imputeacceptrate(b1,0)/=b0;
}
else{
imputeacceptrate(b1,0)*=(b0-1);
imputeacceptrate(b1,0)/=b0;
}

}
}
}


// sample infection time
for (b1=data1.nrow()-1;b1>=0;--b1){
if ((data1(b1,14)==-1)&&(data1(b1,10)==1)){	
accept_pro=1;
IntegerMatrix data11_pro(clone(data11));	
// sample infection time
data11_pro(b1,14)=data1(b1,12)+rand()%(data1(b1,13)+1-data1(b1,12));
newloglik(0)=prior_loglik(p_para(b0,_));
if (newloglik(0)> -9999999){
loglik1pro=loglik(data11_pro,contactmatrix,data2,p_para(b0,_),p_para2(b0,_),SI_vector,inc_vector);
newloglik(0)+=sum(loglik1pro);
loglikeratio=newloglik(0)-temploglik(0);
accept_pro=pow(exp(1),loglikeratio);
}
if(gen_binom(accept_pro)){
loglik1=clone(loglik1pro);		
temploglik(0)=newloglik(0);
data11=clone(data11_pro);
imputeacceptrate(b1,2)*=(b0-1);
imputeacceptrate(b1,2)+=1;
imputeacceptrate(b1,2)/=b0;
}
else{
imputeacceptrate(b1,2)*=(b0-1);
imputeacceptrate(b1,2)/=b0;
}
}
}

// sample symptom time
for (b1=data1.nrow()-1;b1>=0;--b1){
if ((data1(b1,11)==-1)&&(data1(b1,10)==1)){	
accept_pro=1;
IntegerMatrix data11_pro(clone(data11));	
// sample infection time
int bound=data1(b1,13)+15;
if (bound>150){
bound=150;
}
data11_pro(b1,11)=data1(b1,12)+rand()%(bound+1-data1(b1,12));
newloglik(0)=prior_loglik(p_para(b0,_));
if (newloglik(0)> -9999999){
loglik1pro=loglik(data11_pro,contactmatrix,data2,p_para(b0,_),p_para2(b0,_),SI_vector,inc_vector);
newloglik(0)+=sum(loglik1pro);
loglikeratio=newloglik(0)-temploglik(0);
accept_pro=pow(exp(1),loglikeratio);
}
if(gen_binom(accept_pro)){
loglik1=clone(loglik1pro);		
temploglik(0)=newloglik(0);
data11=clone(data11_pro);
imputeacceptrate(b1,3)*=(b0-1);
imputeacceptrate(b1,3)+=1;
imputeacceptrate(b1,3)/=b0;
}
else{
imputeacceptrate(b1,3)*=(b0-1);
imputeacceptrate(b1,3)/=b0;
}
}
}

// sample random effect 
for (b1=data1.nrow()-1;b1>=0;--b1){
if (contactmatrix(b1,0)==-1){	
accept_pro=1;
NumericMatrix data2_pro(clone(data2));	
// sample re
data2_pro(b1,0)=data2(b1,0)+rnorm(0.0,p_para(b0,17));
newloglik(0)=prior_loglik(p_para(b0,_));
if (newloglik(0)> -9999999){
loglik1pro=loglik(data11,contactmatrix,data2_pro,p_para(b0,_),p_para2(b0,_),SI_vector,inc_vector);
newloglik(0)+=sum(loglik1pro);
loglikeratio=newloglik(0)-temploglik(0);
accept_pro=pow(exp(1),loglikeratio);
}
if(gen_binom(accept_pro)){
loglik1=clone(loglik1pro);		
temploglik(0)=newloglik(0);
data2=clone(data2_pro);
imputeacceptrate(b1,4)*=(b0-1);
imputeacceptrate(b1,4)+=1;
imputeacceptrate(b1,4)/=b0;
}
else{
imputeacceptrate(b1,4)*=(b0-1);
imputeacceptrate(b1,4)/=b0;
}
}
}



LL1(b0,1)=temploglik(0)-prior_loglik(p_para(b0,_));

int infindex=0;
for (b1=0;b1<=data11.nrow()-1;++b1){
if (data11(b1,10)==1){	
infection_record(b0,infindex)=data11(b1,14);
++infindex;	
}
}

// record rm
for (b1=data1.nrow()-1;b1>=0;--b1){
if (contactmatrix(b1,0)==-1){
rmrecord(b0,b1)=data2(b1,0);	
}	
}


if (b0%100==1){
//if (b0%1==0){
Rcout << "b0: " << b0 << std::endl;
}

}


return List::create(_[""]=p_para,
_[""]=p_para2,
_[""]=LL1,
_[""]=data11,
_[""]=data2,
_[""]=rmrecord,
_[""]=baselinecount,
_[""]=imputeacceptrate,
_[""]=infection_record);
} 

