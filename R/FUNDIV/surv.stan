data {
    int<lower=0> N; //no. of trees
    real<lower=0> dbh[N]; // diameter in mm
    int<lower=0, upper=1> surv[N]; // survival of each tree, 1 for alive, 0 for dead
    real time[N]; // time interval between two censuses
    real thresh;
    real lowr2;
    real upr2;
}

parameters {
    real<lower=0.01, upper=0.995> K;
    real<lower=-1*thresh, upper=thresh*0.75> p1;
    real<lower=0.0001, upper=0.25> r1;
    real<lower=max(dbh)*1.2, upper=2*(max(dbh))> p2;
    real<lower=lowr2, upper=upr2> r2;
}

model {
    real p;
    real theta;
    for(i in 1:N)
    {
      if(dbh[i] < thresh)
      {
        p = (K) / (1 + exp(-r1 * (dbh[i] - p1)));
        theta = p^time[i];
      }
      else
      {
        p = (K) / (1 + exp(-r2 * (dbh[i] - p2)));
        theta = p^time[i];
      }
    surv[i] ~ bernoulli(theta);
   }
}
