//Feel free to report errors & suggestions: ehsan@stat.ubc.ca
//Last updated: 16-July-2012
//Karim, M. E.; Gustafson, P.; Petkau, J. "Generating survival data for fitting marginal structural Cox models using Stata", Stata Conference, San Diego, July 26–27, 2012
//https://ideas.repec.org/p/boc/scon12/16.html

//clear
//clear matrix
//clear mata


mata

function msm(newx, subjects, tpoints){
rseed(newx)

  true_psi = 0.3 
  L_gamma0 = log(3/7) 
  L_gamma1 = 2 
  L_gamma2 = log(0.5) 
  L_gamma3 = log(1.5)
  A_alpha0 = log(2/7)  
  A_alpha1 = 0.5  
  A_alpha2 = 0.5 
  A_alpha3 = log(4)                      
  constant_point = 30
  scale = 0.01 
  shape = 1 

  output = J(1, 19, .)
  for (id=1; id<=subjects; id++)
  {
    u = uniform(1,1)
    ran_exp = log(1-u)/(-1)
    T0 = (ran_exp/(scale^shape))^(1/shape)
    if (T0 < constant_point) {
                   IT0=1
                } 
                else 
                {
                   IT0=0
                }
    max_T = tpoints
    time_pointer1 = J(1,1,.)
    time_pointer2 = J(1,1,.)

    L = J(1,tpoints,0)
    A = J(1,tpoints,0)
    Y = J(1,tpoints,0)
    Y_prev = J(1,tpoints,0)
    pr_A = J(1,tpoints,0)
    pr_L = J(1,tpoints,0)
    int_factor = J(1,tpoints,0)
    failure_time = J(1,tpoints,0)
    Cense = J(1,tpoints,0)
    Cens  = J(1,tpoints,0)
    pcens = J(1,tpoints,0)
    pcense = J(1,tpoints,0)

    j = 1
    Y_prev[1]=0

    logitL = L_gamma0 + L_gamma1 * IT0
    pr_L[1] = 1/(1+(1/exp(logitL)))
    if (pr_L[1]==1) {
                     L[1] = 1
                  } 
                  else 
                  {
                     successes = uniform(1,1) :< pr_L[1] 
                     L[1] = rowsum(successes) 
                  }
    logitA = A_alpha0 + A_alpha1 * L[1]
    pr_A[1] = 1/(1+(1/exp(logitA)))
    successes = uniform(1,1) :< pr_A[1]
    A[1] = rowsum(successes) 

    int_factor[1] = exp(true_psi*A[1])
    if (T0 > int_factor[1]) {
                        Y[1]=0
                        failure_time[1]=J(1,1,.)
                     } 
                      else 
                     {
                        Y[1]=1
                        failure_time[1] = T0*exp(-true_psi*A[1])
                     }
    j = 1
    Cense[1] = 0
    Cens[1] = 0
    pcens[1] = 0

    for (j=2; j<=tpoints; j++){
      
      if (time_pointer2 != J(1,1,.) ) {
                                  if (j == time_pointer2) {
                                                      j = tpoints
                                                   }
                                }
      if (time_pointer1 != J(1,1,.) ) {
                                  if (j == time_pointer1) {
                                                      j = tpoints
                                                   }
                                }
      Y_prev[j]=Y[j-1]

      logitL = L_gamma0 + L_gamma1 * IT0 + L_gamma2 * A[j-1] + L_gamma3 * L[j-1]
      pr_L[j]=exp(logitL)/(1+exp(logitL))

      if (pr_L[j] == 1) {
                         L[j] = 1
                      } 
                      else {
                              successes = uniform(1,1) :< pr_L[j] 
                              L[j] = rowsum(successes) 
                            }

      logitA = A_alpha0 + A_alpha1 * L[j] + A_alpha2 * L[j-1] + A_alpha3 * A[j-1]
      pr_A[j]=exp(logitA)/(1+exp(logitA))
      successes = uniform(1,1) :< pr_A[j]
      A[j] = rowsum(successes) 

      if (Y_prev[j]==1) {
        A[j]=0
        L[j]=0
                    }


      int_factor[j] = int_factor[j-1] + exp(true_psi*A[j])

      if (T0 > int_factor[j]) {
        Y[j]=0
        failure_time[j]=J(1,1,.)
                        } 
                        else 
                        {
                           Y[j]=1
                           if (Y_prev[j]==1) {
                                            failure_time[j]=failure_time[j-1]
                                          } 
                                          else 
                                          {
                                             failure_time[j]=(j-1)+((T0-int_factor[j-1])*exp(-true_psi*A[j]))
                                           }
                         }


      Cense[j] = Cense[1]
      Cens[j] = 0
      pcens[j] = 0
      pcense[j] = pcense[1] 
    }

    if (failure_time[tpoints] == J(1,1,.)) {
                                    time_pointer1 = tpoints
                                    time_pointer2 = tpoints
                                  } 
                                  else 
                                  {
                                    time_pointer1 = ceil(failure_time[tpoints])
                                    time_pointer2 = ceil(failure_time[tpoints])
                                   }

    active_time = tpoints
      if (time_pointer1 < tpoints){
                        active_time = time_pointer1
                       }  
      if (time_pointer2 < tpoints){
                        active_time = time_pointer2
                       }  
      if (time_pointer2 < time_pointer1){
                            active_time = time_pointer2
                          }
                          else {
                                  active_time = time_pointer1
							   }
					   

    idx = J(1,active_time,0)
    Lx = J(1,active_time,0)
    Ax = J(1,active_time,0)
    Yx = J(1,active_time,0)
    Y_prevx = J(1,active_time,0)
    L_prev = J(1,active_time,0)
    A_prev = J(1,active_time,0)
    A_prevL = J(1,active_time,0)
    pr_A_tx = J(1,active_time,0)
    pr_Lx = J(1,active_time,0)
    int_factorx = J(1,active_time,0)
    Tx = J(1,active_time,0)
    tpointx = J(1,active_time,0)
    tpoint2x = J(1,active_time,0)
    T0x = J(1,active_time,0)
    IT0x = J(1,active_time,0)
    max_Tx = J(1,active_time,0)
    censor1 = J(1,active_time,0)
    censor2 = J(1,active_time,0)
    psi= J(1,active_time,0)
    newxx= J(1,active_time,0)

    for (tpoint=1; tpoint<=active_time; tpoint++)
	  {
      tpoint2 = tpoint-1
      T0x[tpoint]=T0
      IT0x[tpoint]=IT0
      tpointx[tpoint]=tpoint
      tpoint2x[tpoint]=tpoint2
      idx[tpoint]=id
      Lx[tpoint]=L[tpoint]
      int_factorx[tpoint]=int_factor[tpoint]
      pr_Lx[tpoint] = pr_L[tpoint]
      Ax[tpoint]=A[tpoint]
      pr_A_tx[tpoint]=pr_A[tpoint]
      if (tpoint==1) {
                        A_prev[tpoint]=0
                     } 
                     else {
                             A_prev[tpoint]=A[tpoint-1]
                           }
      if (tpoint==1) {
                        L_prev[tpoint]=0
                     } 
                     else {
                             L_prev[tpoint]=L[tpoint-1]
                           }
      A_prevL[tpoint] = A_prev[tpoint]*Lx[tpoint]
                         
      censor1[tpoint] = Cens[tpoint]
      censor2[tpoint] = Cense[tpoint]
                               
      Y_prevx[tpoint] = Y_prev[tpoint]
      psi[tpoint] = true_psi
      newxx[tpoint] = newx

      if (active_time < tpoints) {
        if (tpoint < active_time) {
                                    Yx[tpoint] = Y[tpoint]
                                 } 
                                 else 
                                 {
                                    Yx[tpoint] = 1
                                 }
                                } 
                                else 
                                { 
                                   if (active_time == tpoints) {
                                                                Yx[tpoint] = Y[tpoint]
                                                              } 
                                 }

      Tx[tpoint]=failure_time[max_T]
      max_Tx[tpoint] = max_T

    }


    output = output\(idx', tpointx', tpoint2x', T0x', IT0x', int_factorx', Yx', Y_prevx', Ax', A_prev', Lx', L_prev', A_prevL', pr_A_tx', Tx', max_Tx', pr_Lx', psi', newxx')
    
} 

idx = rows(output)
outputx = output[2..idx,.]
st_matrix("outputx",outputx)
return(outputx)
}

end

