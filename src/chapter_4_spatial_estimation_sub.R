
# required libraries
# require(RcppArmadillo)
if (require(pacman) == FALSE) {
  install.packages("pacman")
}

pacman::p_load(data.table, gridExtra, lattice,
               nplr, tidyverse, raster, rasterVis,
               Rcpp, RcppEigen)


sourceCpp(code = '
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

#include <omp.h>
#include <RcppEigen.h>

// [[Rcpp::export]]
SEXP eigenMapMatMult2(const Eigen::Map<Eigen::MatrixXd> A,
                      Eigen::Map<Eigen::MatrixXd> B,
                      int n_cores){
  Eigen::setNbThreads(n_cores);
  Eigen::MatrixXd C = A * B;
  return Rcpp::wrap(C);
}
')


# cpp functions
cppFunction('
Rcpp::DataFrame combi2inds(const Rcpp::CharacterVector inputVector){
const int len = inputVector.size();
const int retLen = len * (len-1) / 2;
Rcpp::IntegerVector outputVector1(retLen);
Rcpp::IntegerVector outputVector2(retLen);
int indexSkip;
for (int i = 0; i < len; ++i){
    indexSkip = len * i - ((i+1) * i)/2;
    for (int j = 0; j < len-1-i; ++j){
        outputVector1(indexSkip+j) = i+1;
        outputVector2(indexSkip+j) = i+j+1+1;
        }
    }
return(Rcpp::DataFrame::create(Rcpp::Named("xid") = outputVector1,
                          Rcpp::Named("yid") = outputVector2));
};
')



cppFunction('
NumericMatrix crossdist(NumericMatrix m1, NumericMatrix m2) {
  int nrow1 = m1.nrow();
  int nrow2 = m2.nrow();
  int ncol = m1.ncol();
 
  if (ncol != m2.ncol()) {
    throw std::runtime_error("Incompatible number of dimensions");
  }
 
  NumericMatrix out(nrow1, nrow2);
 
  for (int r1 = 0; r1 < nrow1; r1++) {
    for (int r2 = 0; r2 < nrow2; r2++) {
      double total = 0;
      for (int c12 = 0; c12 < ncol; c12++) {
        total += pow(m1(r1, c12) - m2(r2, c12), 2);
      }
      out(r1, r2) = sqrt(total);
    }
  }
 
  return out;
}')



cut_to_lagd <- function(x1){
  lapply(strsplit(gsub("\\(|\\]", "", as.character(x1)), ","),
         function(x2) mean(as.numeric(x2))) %>% unlist
} 


comb_2 <- function(x){
  
  ind <- combi2inds(x)
  return(setDT(list(x[ind$xid], x[ind$yid])))
  
}

comb_diff <- function(x){
  
  ind <- combi2inds(x)
  return(setDT(list(x[ind$xid] - x[ind$yid])))
  
}



# projection CRS
utm_k <- "+proj=tmerc +lat_0=38 +lon_0=127.50289 +k=0.9996 +x_0=1000000 +y_0=2000000 +ellps=bessel +units=m +no_defs"

lambert_aea <- "+proj=laea +lat_0=38 +lon_0=127.5 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

lambert_aea_new <- "+proj=laea +lat_0=34.53333 +lon_0=137.0698 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"


# min-max normalize function
scale_minmax <- function(x) {
  
  (x-min(x))/(max(x)-min(x))
  
}

# min-max normalize function: fixed scale
fixed_scale_minmax <- function(x, fixed_min, fixed_max){
  
  (x - fixed_min)/(fixed_max - fixed_min)
  
}


# lag distance overlapping function
lag_dist <- function(d, buffer, m_factor){
  
  starts <- seq(0, max(d), by = buffer/m_factor)
  ends <- dplyr::lead(starts, m_factor)
  
  lag_interval <- data.frame(lag_st = starts, lag_ed = ends)[1:length(which(!is.na(ends))), ]
  return(lag_interval)
  
}



#############################################################################################
# variogram related functions

## experimental variogram function for spatio-temporal kriging

vge <- function(xi, s_var, val,
                max_range, lag_buffer, overlap_factor,
                valid_n = 30,
                plot = TRUE){
  
  sample_var <- var(xi[, val])
  pair_val <- comb_2(xi[, val])
  
  pair_val[, dist_s := as.numeric(dist(xi[, s_var]))] # add variables[dist_s] on the exist table
  lag_d_s <- lag_dist(pair_val$dist_s, buffer = lag_buffer, m_factor = overlap_factor)
  pair_val[, cut_s :=
             cut(dist_s,
                 c(unique(unlist(lag_d_s)),
                   pair_val[, ceiling(max(dist_s)+lag_buffer)]),
                 right = F)]
  
  pair_val[, n := .(n = length(V1)),
           keyby = .(cut_s)]
  
  # invalid(= n group <= 30) pair filtering
  invalid_pair <- which(pair_val$n <= valid_n)
  
  # experimental variogram calculation
  vg_e = pair_val[-invalid_pair, .(
    cov = cov(V1, V2, use = "everything"),
    # conf_int_cor = cor.test(V1, V2)$conf.int,
    sd_x = sd(V1, na.rm = T),
    sd_y = sd(V2, na.rm = T)),
    keyby = .(cut_s)]
  vg_e[, gamma := sample_var - cov]
  vg_e[, sd_xy := sd_x * sd_y]
  vg_e[, conf_int_cov := sd_x * sd_y]
  vg_e[, lag_s := mean(
    as.numeric(
      gsub("\\[|\\)", "", unlist(str_split(cut_s, ",")))
    )
  ),
       by = .(cut_s)]
  
  # 2 row confidence interval to 1 row
  vg_e[, .(conf_cov_l = conf_int_cov[2], conf_cov_u = conf_int_cov[1]), 
       by = .(cut_s, cov, sd_x, sd_y, sd_xy, lag_s)]
  
  # confirm empirical variogram
  vg_e_lim <- vg_e %>%
    filter(lag_s <= max_range * max(lag_s))
  
  # plotting option
  if(plot == TRUE){
    
    print(
      plot(vg_e_lim$lag_s, vg_e_lim$gamma, pch = 16)
    )
    
    
  }
  
  return(vg_e_lim)
  
  
}




# anisotropic experimental variogram function for xyzt
vge_aniso <- function(xi, s_var, val, aniso = NULL,
                      max_range, lag_buffer, overlap_factor,
                      valid_n = 30,
                      plot = TRUE){
  
  
  aniso <- aniso
  aniso_dist <- paste0("dist_", aniso)
  aniso_cut <- paste0("cut_", aniso)
  aniso_lag <- paste0("lag_", aniso)
  aniso_conf <- paste0(c("conf_cov_l_", "conf_cov_u_"),
                       rep(aniso, each = 2))
  
  sample_var <- var(xi[, val])
  pair_val <- comb_2(xi[, val])
  
  # add variables[dist_x, y, z] on the exist table
  pair_val[, (aniso_dist) := 
             lapply(aniso, function(x) as.numeric(dist(xi[, x])))]
  
  # lag_d_s <- lag_dist(pair_val$dist_x, buffer = lag_buffer, m_factor = overlap_factor)
  
  lag_d_list <- lapply(aniso_dist,
                       function(x){
                         
                         lag_dist(pair_val[[x]],
                                  buffer = lag_buffer,
                                  m_factor = overlap_factor)
                         
                       })
  
  names(lag_d_list) <- aniso_dist
  
  pair_val[, (aniso_cut) :=
             lapply(aniso_dist,
                    function(x){
                      cut(pair_val[[x]],
                          c(unique(unlist(lag_d_list[[x]])),
                            ceiling(max(pair_val[[x]])) + lag_buffer),
                          right = F)
                    })
  ]
  
  
  pair_val[, n := .(n = length(V1)),
           keyby = aniso_cut]
  
  
  # invalid(= n group <= 30) pair filtering
  invalid_pair <- which(pair_val$n <= 10)
  
  # experimental variogram calculation
  vg_e = pair_val[-invalid_pair, .(
    cov = cov(V1, V2, use = "everything"),
    # conf_int_cor = cor.test(V1, V2)$conf.int,
    sd_x = sd(V1, na.rm = T),
    sd_y = sd(V2, na.rm = T)),
    keyby = aniso_cut]
  vg_e[, gamma := sample_var - cov]
  vg_e[, conf_int_gamma := sample_var - sd_x * sd_y]
  
  
  vg_e[, (aniso_lag) := 
         lapply(.SD, function(x){
           mean(
             as.numeric(
               gsub("\\[|\\)", "", unlist(str_split(x, ",")))
             )
           )
         }),
       .SDcols = aniso_cut,
       keyby = aniso_cut]
  
  
  # covariance confidence interval(Not Recommended)
  agg_cols <- c(aniso_cut, "cov", "sd_x", "sd_y", aniso_lag)
  vg_e[, (aniso_conf) := 
         lapply(.SD, function(x){
           c(x[1], x[2])
         }), .SDcols = "conf_int_gamma",
       keyby = agg_cols]
  
  
  # confirm empirical variogram
  vg_e_lim <- vg_e %>%
    filter(if_all(starts_with("lag_"), ~ . <= max_range * max(vg_e$lag_x)))
  
  
  # plotting option
  if(plot == TRUE){
    
    op <- par(no.readonly = T)
    par(mfrow = c(2, 2))
    
    for(pp in aniso_lag){
      
      plot(vg_e_lim[[pp]], vg_e_lim$gamma, pch = 16,
           main = paste0("Gamma For ", pp))
      
    }
    
    par(op)
    
  }
  
  
  return(vg_e_lim)
  
  
}




## theoretical variogram fitting function
vgt <- function(ex_vgm, vgm_list, par0, lower, upper, plot = TRUE){
  
  match_funs <- lapply(vgm_list, match.fun)
  fit <- optim(par0, match_funs[['vgm_s']],
               d = ex_vgm$lag_s, gamma = ex_vgm$gamma,
               return = "error",
               gr=NULL, method = "L-BFGS-B", lower = lower, upper = upper)
  
  new_s <- seq(0, max(ex_vgm$lag_s), length.out = 100)
  
  fitted_gamma <- match_funs[['vgm_s']](d = new_s, gamma = NULL,
                                        pars = fit$par,
                                        return = "value")
  fitted_vgm <- data.frame(new_s, fitted_gamma)
  
  if(plot == TRUE){
    
    # plotting fitted variogram
    # compare with empirical variogram
    plot(ex_vgm$lag_s, ex_vgm$gamma, pch = 16)
    lines(fitted_vgm$new_s, fitted_vgm$fitted_gamma)
    
  }
  
  res <- list(par = fit$par, fitted_vgm = fitted_vgm)
  
  return(res)
  
  
}



#############################################################################################



#############################################################################################
# idw

idw <- function(xi, x0, coord, val,
                order = 2,
                step_size = 100000, n_cores = 8,
                progress = T, silence = F){
  
  # step_size <- step_size
  step_size <- ifelse(step_size > nrow(x0),
                      nrow(x0),
                      step_size)
  max_val = floor(nrow(x0)/step_size)
  
  if(progress == T){
    pb <- txtProgressBar(min = 0, max = max_val, style = 3)
  }
  
  
  res <- rbindlist(
    
    lapply(1:max_val, 
           function(kk){
             
             if(kk != max_val){
               sample_list = (1+(kk - 1) * step_size):(kk * step_size)
             }else{
               sample_list = (1+(kk - 1) * step_size):nrow(x0)
             }
             
             weights <- 1/crossdist(as.matrix(x0[sample_list, coord]),
                                    as.matrix(xi[, coord]))^order
             
             x0_hat <- eigenMapMatMult2(weights,
                                        as.matrix(xi[,val]),
                                        n_cores = n_cores)/rowSums(weights)
             
             if(progress == T){
               setTxtProgressBar(pb, kk)
             }
             
             return(
               
               data.table(est_val = as.numeric(x0_hat))
               
             )
             
             
           }
    )
    
  )
  
}


## idw cv function
idw_cv <- function(xi, coord, val,
                   method = "nfold", n = 10,
                   order = 2,
                   step_size = 100000, n_cores = 8,
                   progress = T, silence = F){
  
  xi <- xi
  
  # shuffled sample index
  fold_id <- cut(sample(1:nrow(xi)), breaks = n, labels = F)
  
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  
  fold_res <- rbindlist(
    
    lapply(1:n,
           
           function(ff){
             
             target_coord <- xi[which(fold_id == ff), coord]
             temp_res1 <- idw(xi = xi[-which(fold_id == ff), ],
                              x0 = target_coord,
                              coord = coord, val = val, order = 2,
                              step_size = step_size, n_cores = n_cores,
                              progress = F, silence = F)
             org_val <- xi[which(fold_id == ff), val]
             temp_res2 <- data.table(target_coord, temp_res1, org_val)
             
             gc(reset = T)
             
             setTxtProgressBar(pb, ff)
             
             return(temp_res2)
             
           }
           
    )
    
  )
  
  return(fold_res)
  
}


#############################################################################################



#############################################################################################

# general kriging core
krige <- function(xi, x0, coord, val, vgm_model = "Sph", vgm_options,
                  step_size = 10000, n_cores = 8,
                  progress = T, silence = F){
  
  obs_h <- crossdist(as.matrix(xi[, coord]), as.matrix(xi[, coord]))
  vgm <- match.fun(vgm_model)
  xs <- as.matrix(rbind(-as.matrix(xi[, coord]), 0))
  xs_t <- t(-xs)
  mat_zero <- matrix(0, nrow = nrow(xs_t), ncol = ncol(xs))
  
  if(silence == F){
    message("\npreparing LHS matrix including matrix inverse...\n please wait...")
  }
  
  mat_L_inv <- vgm(d = obs_h, gamma = NULL,
                   pars = vgm_options$pars,
                   return = "value") %>%
    matrix(., nrow = dim(obs_h)[1], ncol = dim(obs_h)[2]) %T>%
    {assign(x = "ss", value = max(., na.rm = T), envir = .GlobalEnv)} %>%
    {ss - .} %>%
    `diag<-`(., ss) %>%
    cbind(., rep(-1, dim(obs_h)[1])) %>%
    rbind(c(rep(1, dim(obs_h)[2]), 0)) %>%
    cbind(., xs) %>%
    rbind(., cbind(xs_t, mat_zero)) %>%
    solve
  
  # step_size <- step_size
  step_size <- ifelse(step_size > nrow(x0),
                      nrow(x0),
                      step_size)
  max_val = floor(nrow(x0)/step_size)
  
  if(progress == T){
    pb <- txtProgressBar(min = 0, max = max_val, style = 3)
  }
  
  
  res <- rbindlist(
    
    lapply(1:max_val, 
           function(kk){
             
             if(kk != max_val){
               sample_list = (1+(kk - 1) * step_size):(kk * step_size)
             }else{
               sample_list = (1+(kk - 1) * step_size):nrow(x0)
             }
             
             
             target_h <- crossdist(as.matrix(x0[sample_list, coord]),
                                   as.matrix(xi[, coord]))
             
             mat_t <- vgm(d = target_h, gamma = NULL,
                          pars = vgm_options$pars,
                          return = "value") %>%
               matrix(., nrow = dim(target_h)[1], ncol = dim(target_h)[2]) %>%
               {ss - .} %>%
               cbind(., 1) %>%
               cbind(., x0[sample_list,]) %>%
               t
             
             
             coefs <- eigenMapMatMult2(mat_L_inv , mat_t, n_cores = n_cores)
             gc(reset = T)
             
             
             lambda_id <- 1:nrow(xi)
             omega_id <- (nrow(xi)+1):dim(coefs)[1]
             lambda <- coefs[lambda_id,]
             omega <- coefs[omega_id,]
             
             
             est_var1 = diag(eigenMapMatMult2(t(coefs[lambda_id,]),
                                              mat_t[lambda_id,],
                                              n_cores = n_cores))
             est_var2 = diag(eigenMapMatMult2(t(coefs[omega_id,]),
                                              mat_t[omega_id,],
                                              n_cores = n_cores))
             
             if(progress == T){
               setTxtProgressBar(pb, kk)
             }
             
             return(
               
               data.table(est_val = as.numeric(xi[, val] %*% lambda),
                          est_var =  ss - est_var1 + est_var2 )
               
             )
             
             
           }
    )
    
  )
  
  
}

#############################################################################################



#############################################################################################
# cross-validation
krige_cv <- function(xi, coord, val,
                     method = "nfold", n = 10,
                     vgm_model = "Sph", vgm_options,
                     step_size = 10000, n_cores = 8,
                     progress = F, silence = T){
  
  xi <- xi
  
  # shuffled sample index
  fold_id <- cut(sample(1:nrow(xi)), breaks = 10, labels = F)
  
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  
  fold_res <- rbindlist(
    
    lapply(1:n,
           
           function(ff){
             
             target_coord <- xi[which(fold_id == ff), coord]
             temp_res1 <- krige(xi = xi[-which(fold_id == ff), ],
                                x0 = target_coord,
                                coord = coord,
                                val = val,
                                vgm_model = vgm_model,
                                vgm_options = vgm_options,
                                step_size = step_size, n_cores = n_cores,
                                progress = progress, silence = silence)
             org_val <- xi[which(fold_id == ff), val]
             temp_res2 <- data.table(temp_res1, org_val)
             
             gc(reset = T)
             
             setTxtProgressBar(pb, ff)
             
             return(temp_res2)
             
           }
           
    )
    
  )
  
  return(fold_res)
  
}

#############################################################################################




