

# variogram model for spacial ans temporal covariance
Exp <- function(d, gamma, pars, return = "value"){
  
  v_model <- ifelse(d <= pars[3],
                    pars[1] + pars[2]*(1- exp(-3*(d/pars[3]))),
                    pars[1] + pars[2])
  
  if(return == "value"){
    return(v_model)
  }else if(return == "error"){
    return(sum((gamma - v_model)^2))
  }
  
}



Gau <- function(d, gamma, pars, return = "value"){
  
  v_model <- ifelse(d <= pars[3],
                    pars[1] + pars[2]*(1- exp(-3*(d/pars[3])^2)),
                    pars[1] + pars[2])
  
  if(return == "value"){
    return(v_model)
  }else if(return == "error"){
    return(sum((gamma - v_model)^2))
  }
  
}



Sph <- function(d, gamma, pars, return = "value"){
  
  v_model <- ifelse(d <= pars[3],
                    pars[1] + pars[2]*(1.5*(d/pars[3]) -0.5*(d/pars[3])^3),
                    pars[1] + pars[2]) 
  
  if(return == "value"){
    return(v_model)
  }else if(return == "error"){
    return(sum((gamma - v_model)^2))
  }
  
}


