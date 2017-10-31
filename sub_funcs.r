wholeyear_civiltwilight_from_NAVY <- function(Y, lat_dec, lon_dec) {
  ## this is the url for GPS coors.
  ## http://aa.usno.navy.mil/cgi-bin/aa_rstablew.pl?ID=AA&year=2017&task=2&place=&lon_sign=-1&lon_deg=80&lon_min=25&lat_sign=1&lat_deg=43&lat_min=37&tz=&tz_sign=-1 
#  lat_dec <- field_lat; lon_dec <- field_lon;
  lat_sign <- sign(lat_dec); lon_sign <- sign(lon_dec);
  lat_deg <- floor(abs(lat_dec)); lon_deg <- floor(abs(lon_dec));
  lat_min <- floor((abs(lat_dec) - lat_deg) * 60); lon_min <- floor((abs(lon_dec) - lon_deg) * 60);
  tz <- floor(abs(lon_dec) / 15);
  navy_url_pre <- 'http://aa.usno.navy.mil/cgi-bin/aa_rstablew.pl?ID=AA&year=';
  navy_url <- paste(navy_url_pre, Y, '&task=2&place=&lon_sign=', lon_sign, '&lon_deg=', lon_deg, '&lon_min=', lon_min, '&lat_sign=', lat_sign, '&lat_deg=', lat_deg, '&lat_min=', lat_min, '&tz=', tz, '&tz_sign=', lon_sign,  sep = '');
  NAVY_page <- getURL(navy_url);
  NAVY <- read.fwf(textConnection(NAVY_page), header = F, n = 31, skip = 34, widths = c(4, rep(c(5,6), 12)));
  ## next is temperate when there is anouncement for NAVY website
#  NAVY <- read.fwf(textConnection(NAVY_page), header = F, n = 31, skip = 36, widths = c(4, rep(c(5,6), 12)));
  d1 <- paste(Y, '-1-1', sep = ''); dL <- paste(Y, '-12-31', sep = '');
  Ds <- seq(as.Date(d1), as.Date(dL), by = "days");
  DL <- c();
  for (m in 1:12) {
    m1 <- NAVY[,2 * m]; m2 <- NAVY[,2 * m + 1];
    delta_h <- floor(m2/100) - floor(m1/100);
    delta_m <- (m2 - floor(m2 / 100) * 100 - m1 + floor(m1 / 100) * 100) / 60;
    DL <- append(DL,  as.vector(na.omit(round(delta_h + delta_m, digits = 2))));
  }
#  print(NAVY);
  DL_df <- data.frame(date = Ds, DL = DL)
 return (DL_df);
}
####
DayLength_from_NAVY <- function(Y1, Y2, lat, lon, local_file, days) {
# Y1 <- planting_year; Y2 <- ending_year; lat <- field_lat; lon <- field_lon; days <- DAPs;
  DLs <- wholeyear_civiltwilight_from_NAVY(Y1, lat, lon);
  if (Y2 > Y1) {
    DL_2 <- wholeyear_civiltwilight_from_NAVY(Y2, lat, lon); 
    DLs <- rbind(DLs, DL_2);
  }
  DL_window <- DLs[DLs$date %in% days, ];
  write.table(DL_window, file = local_file, sep = "\t", row.name = F, quote = F)
}

#####
TM_from_NOAA <- function(d1, d2, lat, lon, local_file, sts_ghcn, daps) {
#  d1 <- planting_date; d2 <- ending_date; lat <- field_lat; lon <- field_lon; sts_ghcn <- STS_ghcn; daps <- DAPs;
  ds_cutoff <- 0.9 * length(daps);
  field_sts <- meteo_distance(sts_ghcn, lat = lat, lon = lon, radius = 50, limit = 6);
#  print (c(d1, lat, lon));
  st_ids <- unique(as.vector(field_sts$id));
  Tmax_df <-  data.frame(date = daps);
  Tmin_df <- data.frame(date = daps);
  TMs <- meteo_tidy_ghcnd(st_ids[1], var = c('TMAX', 'TMIN'), date_min = d1, date_max = d2);
  #### modify here, get rid of NA records before comparison
  TMs_nNA <- subset(TMs, !is.na(TMs$tmax) & !is.na(TMs$tmin));
  if (nrow(TMs_nNA) < ds_cutoff & length(st_ids) > 1) {
   for (st in st_ids) {
    TMs  <- meteo_tidy_ghcnd(st, var = c('TMAX', 'TMIN'), date_min = d1, date_max = d2);
    TMs_F <- TMs %>%
     mutate(tmax=ifelse(tmax==-9999, NA, 32 + 9 * tmax/50))%>%  # convert to degrees F
     mutate(tmin=ifelse(tmin==-9999, NA, 32 + 9 * tmin/50))%>%  # convert to degrees F
     arrange(date);
    colnames(TMs_F)[c(3,4)] <- paste(st, c(1:2),sep= '_'); 
    Tmax_df <- merge(Tmax_df, TMs_F[,c(2,3)], all.x = T);
    Tmin_df <- merge(Tmin_df, TMs_F[,c(2,4)], all.x = T);
   }
    tmax_m <- round(rowMeans(Tmax_df[,-1], na.rm = T), 2);
    tmin_m <- round(rowMeans(Tmin_df[,-1], na.rm = T), 2);
    
   T_mean <- data.frame(date = daps, TMAX = tmax_m, TMIN = tmin_m)  
  } else {
     TMs_F <- TMs %>%
      mutate(tmax=ifelse(tmax==-9999, NA, 32 + 9 * tmax/50))%>%  # convert to degrees F; noaa record is tenths of degrees C
      mutate(tmin=ifelse(tmin==-9999, NA, 32 + 9 * tmin/50))%>%  # convert to degrees F; noaa record is tenths of degrees C
      arrange(date);
    T_mean <- TMs_F[, c(2, 3, 4)];
    colnames(T_mean) <- c('date', 'TMAX', 'TMIN');
  }
  
  write.table(T_mean, file = local_file, sep = "\t", row.name = F, quote = F);
}

Fill_Missning_TM <- function(Tx, env) {
# Tx <- DL_TM$TMAX;
 NA_inds <- which(is.na(Tx));
 A_inds <- which(!is.na(Tx))
 for (NA_ind in NA_inds ) {
   if (NA_ind == 1) { 
     Tx[NA_ind] <- Tx[A_inds[1]]
     } else if (NA_ind >= max(A_inds)) {
       Tx[NA_ind] <- Tx[max(A_inds)]
       } else {
         pre_ind <- A_inds[max(which(A_inds < NA_ind))]; 
         suff_ind <- A_inds[min(which(A_inds > NA_ind))];
         Tx[NA_ind] <- mean(Tx[c(pre_ind, suff_ind)] );
        }
 }
 return(Tx);
}

Adjusting_TM_4_GDD <- function(Tmax, Tmin, threshold, t_base, t_max1, t_max2) {
# Tmax <- Tmax[1:10]; Tmin <- Tmin[1:10]; threshold <-  Haun_threshold
 Tmax[Tmax < t_base] <- t_base; Tmin[Tmin < t_base] <- t_base;
 if (Tmax[1] > t_max2) {Tmax[1] <- t_max2}; 
 if (Tmin[1] > t_max2) {Tmin[1] <- t_max2};
 if (threshold > 0) {
   gdd_cum <- (Tmax[1] + Tmin[1]) / 2 - t_base;
   for (i in 2:length(Tmax)) {
    if (gdd_cum < threshold) { t_max0 <- t_max2} else {t_max0 <- t_max1};
    if (Tmax[i] > t_max0) { Tmax[i] <- t_max0 };
    if (Tmin[i] > t_max0) { Tmin[i] <- t_max0 };
    gdd_cum <- gdd_cum + (Tmax[i] + Tmin[i]) / 2 - t_base;
    }    
  } else {
     Tmax[Tmax > t_max1] <- t_max1; Tmin[Tmin > t_max1] <- t_max1;
     }
 gdds <- round((Tmax + Tmin) / 2 - t_base, 4)
 return (gdds)
} 

Compile_PTT_PTR <-  function(exp_dir, env_meta_info, exp_s_year, exp_e_year, searching_daps,ptt_ptr_file) {
 sp_env_dir <- paste(exp_dir, 'envs/', sep = '');      if (!dir.exists(sp_env_dir))  { dir.create(sp_env_dir)};
 sp_isd_dir  <- paste(sp_env_dir, 'isd/', sep = '');  if (!dir.exists(sp_isd_dir))  { dir.create(sp_isd_dir)};
 sp_ghcn_dir <- paste(sp_env_dir, 'ghcn/', sep = ''); if (!dir.exists(sp_ghcn_dir)) { dir.create(sp_ghcn_dir)};
 sp_navy_dir <- paste(sp_env_dir, 'navy/', sep = ''); if (!dir.exists(sp_navy_dir)) { dir.create(sp_navy_dir)};
 env_meta_info <- env_meta_info_0;
 lat_range <- range(env_meta_info$lat, na.rm = T); lon_range <- range(env_meta_info$lon, na.rm = T);
 local_ghcn_st_file <- paste('D:/0GbE/all_ghcn_stations', sep = '');
 local_ghcn_target_st_file <- paste(sp_ghcn_dir, '0target_ghcn_stations', sep = '');
 if (!file.exists(local_ghcn_target_st_file)) {
   if (!file.exists(local_ghcn_st_file)) { ghcn_all_sts <- ghcnd_stations()} else {ghcn_all_sts <- read.csv(local_ghcn_st_file);}
   sts_target <- dplyr::filter(ghcn_all_sts, first_year <= exp_s_year & last_year >= exp_e_year 
                                             & between(latitude, lat_range[1] - 2, lat_range[2] + 2) 
                                             & between(longitude,lon_range[1] - 2, lon_range[2] + 2)
                                             & element %in% c("TMAX", "TMIN")
                                             );
   write.csv(sts_target, local_ghcn_target_st_file  )
 }
 STS_ghcn <- read.csv(local_ghcn_target_st_file);  
 ghcnd_clear_cache()
 PTT_PTR <- matrix(ncol = 8, nrow = 0); 
 for (e_i in 1:nrow(env_meta_info)) { ## 
   Sys.sleep(1);
  env_code <- env_meta_info$env_code[e_i];
  field_lat <- env_meta_info$lat[e_i]; field_lon <- env_meta_info$lon[e_i]; 
  planting_date <- env_meta_info$PlantingDate[e_i];
  planting_year <- year(planting_date);
  DAPs <- seq(as.Date(planting_date), length.out = searching_daps, by = "day");
  ending_date <- DAPs[searching_daps]; 
  ending_year <- year(ending_date);
  local_dl_file <- paste(sp_navy_dir, env_code, '_DL_', searching_daps, 'DAPs', sep = '')
#  local_tm_file <- paste(sp_ghcn_dir, env_code, '_TM_', searching_daps, 'DAPs', sep = '' );
  local_tm_file <- paste(sp_isd_dir, env_code, '_TM_', searching_daps, 'DAPs', sep = '' );
  if (!file.exists(local_dl_file)) {DayLength_from_NAVY(planting_year, ending_year, field_lat, field_lon, local_dl_file, DAPs ) }
  if (!file.exists(local_tm_file)) {TM_from_NOAA       (planting_date, ending_date, field_lat, field_lon, local_tm_file, STS_ghcn, DAPs)} ;
  DL <- read.table(local_dl_file, header = T, sep = "\t", stringsAsFactors = F);
#  DL$DL[which(DL$DL > 16)] <- 16;
  TM_0 <- read.table(local_tm_file, header = T, sep = "\t", stringsAsFactors = F);
 # TM <- TM_0[,-1];  
  DL_TM <- merge(DL, TM_0, all.x = T) ;
   
  DL_TM$TMAX <- Fill_Missning_TM(DL_TM$TMAX, env);
  DL_TM$TMIN <- Fill_Missning_TM(DL_TM$TMIN, env);
  Tmax <- DL_TM$TMAX; Tmin <- DL_TM$TMIN;
  GDDs <- Adjusting_TM_4_GDD(Tmax, Tmin, Haun_threshold, t_base, t_max1, t_max2);
  PTTs <- round(GDDs * DL_TM$DL, 4);
  PTRs <- round((DL_TM$TMAX - DL_TM$TMIN) / DL_TM$DL, 4);
  PTSs <- round(((DL_TM$TMAX^2) - (DL_TM$TMIN^2)) * (DL_TM$DL^2), 4);
  env_codes <- rep(env_code, searching_daps);
  
  t_df <- data.frame(env_code = env_code, date = DAPs, TMAX = DL_TM$TMAX, TMIN = DL_TM$TMIN, DL = DL_TM$DL, GDD = GDDs, PTT = PTTs, PTR = PTRs, PTS = PTSs)
  PTT_PTR <- rbind(PTT_PTR, t_df);
 }
 
 write.table(PTT_PTR, file = ptt_ptr_file, sep = "\t", quote = F, row.name = F)

}
#########
Pairwise_trait_env_distribution_plot <- function(exp_trait, exp_trait_dir, trait, all_env_codes, env_meta_info) {
 env_mean_trait <- aggregate(x = exp_trait$Yobs, by = list(env_code = exp_trait$env_code), mean);
 n_envs <- nrow(env_mean_trait);
 trait_dist_pdf_file <- paste(exp_trait_dir, trait, '_dist_', n_envs, 'envs.pdf', sep = '');
 pairwise_pdf_file <- paste(exp_trait_dir, trait, '_pairwise_dis', n_envs, 'envs.pdf', sep = '');
 mse_file <- paste(exp_trait_dir, n_envs, 'Env_meanY_MSE', sep = '');

 colnames(env_mean_trait)[2] <- 'meanY';
 env_mean_trait <- env_mean_trait[order(env_mean_trait$meanY),];
 n_obs_env <- c(); quantile_1 <- c(); quantile_3 <- c(); key_para <- c(); key_para2 <- c();
 for (k in 1:nrow(env_mean_trait)) {
   env_data <- subset(exp_trait, exp_trait$env_code == env_mean_trait[k,1]);
   quantiles <- quantile(env_data$Yobs, na.rm = T);
   quantile_1[k] <- quantiles[2];
   quantile_3[k] <- quantiles[4];
   n_obs_env[k] <- length(which(!is.na(env_data$Yobs)))
 }
 env_mean_trait$q1 <- quantile_1; env_mean_trait$q3 <- quantile_3; env_mean_trait$n_obs <- n_obs_env;
 
# line_cnt <- aggregate(x = exp_trait$Yobs, by = list(line_code = exp_trait$line_code),  length);
# line_cnt <- line_cnt[order(line_cnt$x,decreasing = T),]
 line_codes <- unique(as.vector(exp_trait$line_code)); 
 
# line_colors <- terrain.colors(length(line_code));# colorRampPalette(brewer.pal(9,"Blues"))(250); #(terrain.colors(250))
 gray_alpha <- rgb(128, 128, 128, alpha = 35, maxColorValue = 255);
 poly_alpha <- rgb(238, 130, 238, alpha = 55.5, maxColorValue = 255);

 env_codes <- as.vector(env_mean_trait$env_code);
 line_by_env_df <- data.frame(line_code = line_codes);
 for (e_i in 1:nrow(env_mean_trait)) {
   e <- env_codes[e_i];
   e_trait <- subset(exp_trait, exp_trait$env_code == e);
#   nonNAs <- length(which(!is.na(e_trait[,3])))
#   if (nonNAs > (0.5 * length(line_codes))) {
    colnames(e_trait)[3] <- e;
    line_by_env_df <- merge(line_by_env_df, e_trait[,c(1,3)], all.x = T)
#   }
 }
# line_by_env_df$cnt <- apply(line_by_env_df[,-1], 1, function(x) length(na.omit(x)) );
# line_by_env_df <- line_by_env_df[order(line_by_env_df$cnt, decreasing = T),]
 write.table(line_by_env_df, file = paste(exp_trait_dir, 'LbE_table', sep = ''), sep = "\t", row.names = F, quote = F);
 lm_residuals <- data.frame(env_code = env_mean_trait$env_code);
  for (l in line_codes) {
    line_trait_0 <- subset(exp_trait, exp_trait$line_code == l & !is.na(exp_trait$Yobs));
    if (nrow(line_trait_0) > 0) {
     line_trait_0 <-  merge(env_mean_trait, line_trait_0[,c(2,3)])
     line_lm <- lm(line_trait_0$Yobs ~ line_trait_0$meanY);
     df1 <- data.frame(env_code = line_trait_0$env_code, line_code = round((line_lm$residuals)^2, 3));
     colnames(df1)[2] <- l;
     lm_residuals <- merge(lm_residuals, df1, all.x = T)
    }
  }
  df2 <- data.frame(env_code = lm_residuals[,1], errors = rowMeans(lm_residuals[,-1], na.rm = T));
  df2 <- merge(df2, env_mean_trait);
  write.table(lm_residuals, file = mse_file, sep = "\t", row.names = F, quote = F);

 pdf(pairwise_pdf_file, width= 3,height= 3,pointsize=6);
 corrgram(as.matrix(line_by_env_df[,-1]), order=TRUE, lower.panel=panel.ellipse, pch = 19, upper.panel=panel.pie);
 dev.off();
 
 pdf(trait_dist_pdf_file, width= 6,height= 2,pointsize=6);

 layout(matrix(c(1:3), 1, 3, byrow = T))
 
 env_geo_order_df <- merge(env_mean_trait, env_meta_info);
 env_geo_order_df <- env_geo_order_df[order(env_geo_order_df$lat, env_geo_order_df$lon, env_geo_order_df$PlantingDate),];
 env_geo_order <- match(env_mean_trait$env_code, env_geo_order_df$env_code);
# env_mean_trait <- merge(env_meta_info, env_mean_trait);
 par(mar = c(5.0, 2.0, 1, 0.5) , mgp = c(1, 0.1, 0), tck = -0.01, cex.axis = .7, cex.lab = .8, family = "mono");
 plot(0, 0, col = "white", xlim = range(env_geo_order), ylim = range(exp_trait$Yobs, na.rm = T),  ylab = trait,  xlab = '', fg = "gray50", xaxt = "n");
 for (i in 1:nrow(line_by_env_df)) {
  df4 <- data.frame(env_code = env_mean_trait$env_code,  env_order = env_geo_order, Yobs = as.numeric(line_by_env_df[i, -1]));
  df4 <- df4[!is.na(df4$Yobs),];
  df4 <- df4[order(df4$env_order),];
  points(df4$env_order, df4$Yobs, col = gray_alpha, type = "l", pch = 19, lwd = .3)
  points(df4$env_order, df4$Yobs, col = gray_alpha,  pch = 19, cex = .3)
  
#   points(env_mean_trait$meanY, as.vector(line_by_env_df[i, -1]), col = gray_alpha, type = "l", pch = 19, lwd = .3)
#   points(env_mean_trait$meanY, as.vector(line_by_env_df[i, -1]), col = gray_alpha,  pch = 19, cex = .3)
 }
  points(c(1:nrow(env_geo_order_df)), env_geo_order_df$meanY, col = env_cols[match(as.vector(env_geo_order_df$env_code), all_env_codes )], cex = .8)
  points(c(1:nrow(env_geo_order_df)), env_geo_order_df$meanY, col = "black", type = "l", lwd = .5)

 mtext(env_geo_order_df$env_code, side = 1, at = c(1:nrow(env_geo_order_df)), las = 2, line = 0.5, cex = .6 )

 x_tick <- diff(env_mean_trait[,2]) / 50;
 par(mar = c(2.0, 2.0, 1, 0.5) , mgp = c(1, 0.1, 0), tck = -0.01, cex.axis = .7, family = "mono");
 plot(0, 0, col = "white", xlim = range(env_mean_trait$meanY), ylim = range(exp_trait$Yobs, na.rm = T), cex.lab = .9,  ylab = trait,  xlab = 'Population mean', fg = "gray50");
 for (i in 1:nrow(line_by_env_df)) {
   df3 <- data.frame(meanY = env_mean_trait$meanY, Yobs = as.numeric(line_by_env_df[i, -1]));
   df3 <- df3[!is.na(df3$Yobs),];
   points(df3$meanY, df3$Yobs, col = gray_alpha, type = "l", pch = 19, lwd = .3)
   points(df3$meanY, df3$Yobs, col = gray_alpha,  pch = 19, cex = .3)
   
#   points(env_mean_trait$meanY, as.vector(line_by_env_df[i, -1]), col = gray_alpha, type = "l", pch = 19, lwd = .3)
#   points(env_mean_trait$meanY, as.vector(line_by_env_df[i, -1]), col = gray_alpha,  pch = 19, cex = .3)
 }
 abline(a = 0, b = 1, lty = 2, col = "grey")
  polygon(c(env_mean_trait[,2], rev(env_mean_trait[,2])), c( env_mean_trait$q1 , rev(env_mean_trait$q3)), col = poly_alpha, border = "NA")
 points(env_mean_trait$meanY, env_mean_trait$meanY, col = "black", cex = .4)
 legend("topleft", as.vector(env_mean_trait$env_code), col = env_cols[match(as.vector(env_mean_trait$env_code), all_env_codes )], pch = 19, bty = "n", cex = .8)
  
 for (k in 1:nrow(env_mean_trait)) {
   env_data <- subset(exp_trait, exp_trait$env_code == env_mean_trait[k,1]);
   boxplot(env_data$Yobs, add = T, boxwex = x_tick * 5 * 2, at = env_mean_trait$meanY[k], cex = .3, border = env_cols[match(as.vector(env_data$env_code), all_env_codes )], lwd = .4, boxlwd = .7, medlwd = .5, yaxt = "n")
 }
 
 par(mar = c(4.0, 2.0, 1, 0.5) , mgp = c(0.75, 0.1, 0), tck = -0.01, cex.axis = .7, cex.lab = .8,  family = "mono");
 plot(df2$meanY, df2$errors, ylab = 'MSE', xlab = '', xaxt = "n", pch = 19, col = env_cols[match(as.vector(df2$env_code), all_env_codes )]);
 mtext(df2$env_code, side = 1, at = df2$meanY, las = 2, line = 0.5, cex = .6 )
# axis(1, at = df2$meanY, lab = df2$env_code, las = 2)

 dev.off()
}

#################
Exhaustive_search <- function(env_mean_trait, env_paras, searching_daps, exp_trait_dir, FTdaps, trait, p, dap_x, dap_y) {
# env_paras <- PTT_PTR; FTdaps <- exp_traits$FTdap; p <- 1; dap_x <- searching_daps; dap_y <- searching_daps;
 pop_cor_file <- paste(exp_trait_dir, trait, '_', nrow(env_mean_trait), 'Envs_PTTPTR_cor', sep = '');
 exs_pdf_file <- paste(exp_trait_dir, 'MaxR_',trait, '_', nrow(env_mean_trait), 'Envs.pdf', sep = ''); 
# env_paras[,9] <- (env_para[,3]^2 - env_paras[,4]^2) / env_paras[,5];

 if (!file.exists(pop_cor_file)) {
  dap_win <- searching_daps * searching_daps  / 2;
      
    pop_cors_matrix <- matrix(ncol = 14, nrow = dap_win * 1);
    colnames(pop_cors_matrix) <- c("pop_code", 'Day_x', 'Day_y', 'window', 'R_DL', 'R_GDD', 'R_PTT', 'R_PTR', 'R_PTS', 'nR_DL', 'nR_GDD', 'nR_PTT', 'nR_PTR', 'nR_PTS');
    n <- 0;
      for (d1 in 1:(dap_y -1)) {
        for (d2 in (d1 + 1):dap_y) {
         days <- c(d1:d2); 
         env_facts_matrix <- matrix(nrow = nrow(env_mean_trait), ncol = 5);
         for (e_i in 1:nrow(env_mean_trait)) {
           e <- env_mean_trait$env_code[e_i];
           env_para <- subset(env_paras, env_paras$env_code == e);
           env_mean <- colMeans(env_para[days, (c(8, 9, 10, 11, 12) - 3)]); ### DL, GDD, PTT, PTR, PTS
           env_facts_matrix[e_i,] <- env_mean;
         }
         n <- n + 1;
         ### leave one environment out and get the median correlation
         Ymean_envPara <- cbind(env_facts_matrix, env_mean_trait$meanY)
         loo_rs_matrix <- matrix(nrow = nrow(Ymean_envPara)+ 0, ncol = 5);
         for (k in 1:5) {
          for (e_x in c(1:nrow(Ymean_envPara))) {
            t_matrix <- Ymean_envPara[-e_x,];
            loo_rs_matrix[e_x, k] <- round(cor(t_matrix[,6], t_matrix[,k]), digits = 4)
           }
          }
          rs <- apply(loo_rs_matrix, 2, median);
         pop_cors_matrix[n, ] <- c(p, d1, d2, d2 - d1, rs, 0 - rs);
        }
      }
    pop_cors_matrix <- pop_cors_matrix[1:n,]
    write.table(pop_cors_matrix, file = pop_cor_file, sep = "\t", row.names = F, quote = F);

 }

 pop_cors <- read.table(pop_cor_file, header = T, sep = "\t");
 pop_cor <- subset(pop_cors, pop_cors$pop_code == p);
 pdf(exs_pdf_file,width= 5,height= 2,pointsize=6)
 layout(matrix(c(1:10), 2, 5, byrow = T))
 
 for (k in 1:10) {
   pop_cor_0 <- subset(pop_cors, pop_cors$pop_code == p); 
   pop_cor <- pop_cor_0[,c(1:4, k + 4)];
   colnames(pop_cor)[5] <- 'R';
   pop_cor <- pop_cor[order(pop_cor$R),];
   
   xs <- pop_cor$Day_x;  ys <-  pop_cor$Day_y;  mid_R <- median(pop_cor$R);
#   pop_cor_L <- subset(pop_cor, pop_cor$R <= mid_R); cor_range_L <- range(pop_cor_L$R);
#   cell_col_L <- floor((pop_cor_L$R - min(pop_cor_L$R)) / diff(cor_range_L) * col_wdw / 2) + 1;
#  
#   pop_cor_G <- subset(pop_cor, pop_cor$R > mid_R); cor_range_G <- range(pop_cor_G$R);
#   cell_col_G <- ceiling((pop_cor_G$R - min(pop_cor_G$R)) / diff(cor_range_G) * col_wdw / 2) + 12;
#   cell_col <- c(cell_col_L, cell_col_G);
#   pop_cor$cell_col <- cell_col; 
    
   cell_col <- floor(pop_cor$R * 12) + 13; ### the same color scale
   pop_cor$cell_col <- cell_col; 
 
#   pop_cor_m <- subset(pop_cor, pop_cor$window > 4 & pop_cor$window < 50 ); ##& pop_cor$Day_x > (window_ref$dap_s - 10) & pop_cor$Day_y < (window_ref$dap_e + 10));
#   max_R <- pop_cor_m[which.max(pop_cor_m$R)[1], ];
   max_R <- pop_cor[which.max(pop_cor$R)[1], ];
   
   par(mar = c(0.5, 1.0, 1, 0.5) , mgp = c(0.05, 0.1, 0), tck = -0.01, bty = "n");
   plot(-50, -50, xlim = c(0, dap_x), ylim = c(0, dap_y), col = "white",  xlab = '', xaxt = "n", yaxt = "n", ylab = 'Days after planting', bty = "n", cex.lab = .4);
   arrows(-1, 10, -1, dap_y - 10, length = 0.05, angle = 15, lwd = .5,  col = "grey59");
   mtext(c(1, 50, 100, dap_y), side = 2, at = c(1,50, 100, dap_y), line = -1, cex = .6)
  
   rect(xs - 0.5, ys - 0.5, xs + 0.5, ys + 0.5, col = col_palette[pop_cor$cell_col], border = "NA")
   rect(max(pop_cor$Day_x) - 0.5, max(pop_cor$Day_y) - 0.5, max(pop_cor$Day_x) + 0.5, max(pop_cor$Day_y) + 0.5, border = "NA", col = "white", lwd = 0.001)
#   legend("bottom", Div_Fnd_lab, bty = "n", cex = .6);
  
   arrows(10, dap_y + 4, dap_x - 10, dap_y + 4, angle = 15, length = 0.05, lwd = .5, col = "grey59")
   mtext("Days after planting", side = 3, at = dap_x / 2, line = -0.1, cex = .4)
   mtext(c(1, 50, 100, dap_y), side = 3, at = c(1, 50, 100, dap_y), line = -1.1, cex = .6)
   arrows(max_R$Day_x + 4,  max_R$Day_y - 4,  max_R$Day_x,  max_R$Day_y, length = 0.05, angle = 15, lwd = .5, col = "grey59")
   
   box_ys <- seq(1, 50, by = 2); box_xs <- rep(dap_x - 15, 25); 
   rect(box_xs - .5 * 2, box_ys - 0.5 * 2, box_xs + 0.5 * 2, box_ys + 0.5 * 2, border = "NA", col = col_palette)
   text(dap_x - 10 - 5, 52, 'r', cex = .5);
   r_lab_top <- 1; r_lab_mid <- 0; r_lab_bottom <- -1; max_r_lab <- paste( 'r = ', sprintf( "%.3f", max_R$R), sep = '');
   if (k >= 5) { r_lab_top <- -1; r_lab_bottom <- 1; max_r_lab <- paste( 'r = ', sprintf( "%.3f", 0 - max_R$R), sep = '');}
   legend(max_R$Day_x - 4 , max_R$Day_y - 4 , c(paste( max_R$Day_x, ' to ', max_R$Day_y, ' DAP', sep = ''), max_r_lab),  cex = .6, bty = "n");
   text(dap_x - 10 + 3, 50, r_lab_top, cex = .5)
   text(dap_x - 10 + 3, 27, r_lab_mid, cex = .5);
   text(dap_x - 10 + 3, 1,  r_lab_bottom, cex = .5)
   boxplot(FTdaps,   at = 145,  add = TRUE, xaxt = "n", yaxt = "n", xlab = '', ylab = '', width = 10, pch = 19, cex = .3, boxwex = 4, lwd = .4, col = "gold", border = "grey");
   boxplot(FTdaps,   at = 1, horizontal = T, add = TRUE, xaxt = "n", yaxt = "n", xlab = '', ylab = '', width = 10, pch = 19, cex = .3, boxwex = 4, lwd = .4, col = "gold", border = "grey");
   text(mean(FTdaps), 5, 'Days to anthesis', cex = .5)
   text(mean(FTdaps), 10, paste('Trait: ', trait, sep = ''), cex = .6)
 }
 dev.off()
 

}

##### 
Exhaustive_search_animation <- function(env_mean_trait, PTT_PTR, searching_daps, exp_trait_dir, FTdaps, trait, p, dap_x, dap_y) {
  ani_pdf_file <- paste(exp_trait_dir, trait, '_', nrow(env_mean_trait), 'ExSch_ani', sep = '')
  pop_cor_file <- paste(exp_trait_dir, trait, '_', nrow(env_mean_trait), 'Envs_PTTPTR_cor', sep = '');
  exs_pdf_file <- paste(exp_trait_dir, 'MaxR_',trait, '_', nrow(env_mean_trait), 'Envs.pdf', sep = ''); 

  if (!file.exists(pop_cor_file)) {
  saveLatex (
  { par(mar = c(1, 1, 1, 0.5), mgp = c(2, 0.5, 0), tcl = -0.3, cex.axis = 0.8, cex.lab = 0.8, cex.main = 1);
   dev.hold()   
   layout(matrix(c(1,2), 1, 2, byrow = TRUE));
   dap_win <- searching_daps * searching_daps  / 2;
   pop_cors_matrix <- matrix(ncol = 5, nrow = dap_win);
   colnames(pop_cors_matrix) <- c("pop_code", 'Day_x', 'Day_y', 'R_PTT', 'R_PTR'); ## pop_code is reserved in case for multiple familes like NAM
   n <- 0;
   for (d1 in 1: (searching_daps - 1)) { ## 
     for (d2 in (d1 + 1):(searching_daps - 1)) { ## searching_daps
      days <- c(d1:d2);  env_ptts <- c(); env_ptrs <- c();
      for (e in env_mean_trait[,1]) {
        env_para <- subset(PTT_PTR, PTT_PTR$env_code == e);
        env_ptt_mean <- mean(env_para$PTT[days]); env_ptts <- append(env_ptts, env_ptt_mean)
        env_ptr_mean <- mean(env_para$PTR[days]); env_ptrs <- append(env_ptrs, env_ptr_mean) 
      }
      n <- n + 1;
      r_ptt <-  round(cor(env_mean_trait$meanY, env_ptts), digits = 4);
      r_ptr <-  round(cor(env_mean_trait$meanY, env_ptrs), digits = 4);
       
      pop_cors_matrix[n, ] <- c(p, d1, d2, r_ptt, r_ptr);
      if (d1 > 25 & d1 < 50 & d2 < 50)  {
       par(mar = c(2.5, 1.5, 0.75, 0.5) , mgp = c(.75, 0.025, 0), tck = -0.005, family = "mono");
       plot(env_mean_trait$meanY, env_ptts, col = env_cols[match(env_mean_trait$env_code, all_env_codes)], xlab = 'pop mean', ylab = 'PTT', pch = 19, cex = .7);
       abline(lm(env_ptts ~ env_mean_trait$meanY), lty = 2, col = "gray");
       legend("top", paste('DAPs: ', d1, ' - ', d2, ' r = ', r_ptt, sep = ''), bty = "n", cex = .7)
       par(mar = c(2.5, 1.5, 0.75, 0.5) , mgp = c(.75, 0.025, 0), tck = -0.005, family = "mono");
       plot(env_mean_trait$meanY, env_ptrs, col = env_cols[match(env_mean_trait$env_code, all_env_codes)], xlab = 'pop mean', ylab = 'PTR', pch = 19, cex = .7);
       abline(lm(env_ptrs ~ env_mean_trait$meanY), lty = 2, col = "gray");
       legend("top", paste('DAPs: ', d1, ' - ', d2, ' r = ', r_ptr, sep = ''), bty = "n", cex = .7)
       ani.pause()
       }
     }
   }
   pop_cors_matrix <- pop_cors_matrix[1:n,]
   write.table(pop_cors_matrix, file = pop_cor_file, sep = "\t", row.name = F, quote = F);
 
  },
  img.name = ani_pdf_file, ani.opts = "controls,loop,width=0.95\\textwidth", outdir = exp_trait_dir,
  latex.filename = ifelse(interactive(), "ExS.tex", ""), ani.width = 6, ani.height = 3, loop = T,
  interval = 0.2, nmax = 1250, ani.dev = "pdf", ani.type = "pdf", 
  documentclass = paste("\\documentclass{article}","\\usepackage[papersize={6in,3.5in},margin=0.2in]{geometry}", sep = "\n")
  )

 }

 pop_cors <- read.table(pop_cor_file, header = T, sep = "\t");
 pop_cor <- subset(pop_cors, pop_cors$pop_code == p);
 ptt_ptr_ind <- c(4, 4, 5, 5); pm_sign <- c(1, -1, 1, -1);
 ptt_ptr_r_lab <- c('PTT', '0 - PTT', 'PTR', '0 - PTR');
# dap_x <- 150; dap_y <- 150;
 pdf(exs_pdf_file,width= 4,height= 4,pointsize=6)

 layout(matrix(c(1:4), 2, 2, byrow = T))
 for (k in 1:4) {
   pop_cor$R <- pm_sign[k] * pop_cor[, ptt_ptr_ind[k]];
   pop_cor <- pop_cor[order(pop_cor$R),];
   xs <- pop_cor$Day_x;  ys <-  pop_cor$Day_y;  mid_R <- median(pop_cor$R);
   pop_cor_L <- subset(pop_cor, pop_cor$R <= mid_R); cor_range_L <- range(pop_cor_L$R);
   cell_col_L <- floor((pop_cor_L$R - min(pop_cor_L$R)) / diff(cor_range_L) * col_wdw / 2) + 1;
  
   pop_cor_G <- subset(pop_cor, pop_cor$R > mid_R); cor_range_G <- range(pop_cor_G$R);
   cell_col_G <- ceiling((pop_cor_G$R - min(pop_cor_G$R)) / diff(cor_range_G) * col_wdw / 2) + 12;
   cell_col <- c(cell_col_L, cell_col_G);
   pop_cor$cell_col <- cell_col; 
    
   maxR_window <- pop_cors
   maxR_window <- subset(pop_cor, pop_cor$Day_x < 100 & pop_cor$Day_y <= 100);
 
   max_R <- maxR_window[which.max(maxR_window $R)[1], ];
#   screen(k)
   par(mar = c(0.5, 1.0, 1, 0.5) , mgp = c(0.05, 0.1, 0), tck = -0.01, bty = "n", family = "mono");
   plot(-50, -50, xlim = c(0, dap_x), ylim = c(0, dap_y), col = "white",  xlab = '', xaxt = "n", yaxt = "n", ylab = 'Days after planting', bty = "n", cex.lab = .6);
   arrows(-1, 10, -1, dap_y - 10, length = 0.05, angle = 15, lwd = .5,  col = "grey59");
   mtext(c(1, 50, 100, dap_y), side = 2, at = c(1,50, 100, dap_y), line = -1, cex = .6)
  
   rect(xs - 0.5, ys - 0.5, xs + 0.5, ys + 0.5, col = col_palette[pop_cor$cell_col], border = "NA")
   rect(max(pop_cor$Day_x) - 0.5, max(pop_cor$Day_y) - 0.5, max(pop_cor$Day_x) + 0.5, max(pop_cor$Day_y) + 0.5, border = "NA", col = "white", lwd = 0.001)
 #  legend("bottom", Div_Fnd_lab, bty = "n", cex = .6);
  
   arrows(10, dap_y + 4, dap_x - 10, dap_y + 4, angle = 15, length = 0.05, lwd = .5, col = "grey59")
   mtext("Days after planting", side = 3, at = dap_x / 2, line = -0.1, cex = .6)
    mtext(c(1, 50, 100, dap_y), side = 3, at = c(1, 50, 100, dap_y), line = -1.1, cex = .6)
   arrows(max_R$Day_x + 4,  max_R$Day_y - 4,  max_R$Day_x,  max_R$Day_y, length = 0.05, angle = 15, lwd = .5, col = "grey59")
   legend(max_R$Day_x - 4 , max_R$Day_y - 4 , c(paste( max_R$Day_x, ' to ', max_R$Day_y, ' DAP', sep = ''), paste( 'r = ', sprintf( "%.3f", pm_sign[k] * max_R$R), sep = '')),  cex = .6, bty = "n")
   
   box_ys <- seq(1, 50, by = 2); box_xs <- rep(dap_x - 15, 25); 
   rect(box_xs - .5 * 2, box_ys - 0.5 * 2, box_xs + 0.5 * 2, box_ys + 0.5 * 2, border = "NA", col = col_palette)
   text(dap_x - 10 - 10, 52 + 3, paste('r: ', trait, ' & ', ptt_ptr_r_lab[k], sep = ''), cex = .5);
   text(dap_x - 10 + 3, 50, round(pm_sign[k] * max(pop_cor$R), 2), cex = .5)
   text(dap_x - 10 + 3, 27, round(pm_sign[k] * mid_R, 2), cex = .5);
   text(dap_x - 10 + 3, 1,  round(pm_sign[k] * min(pop_cor$R), 2), cex = .5)
   boxplot(FTdaps,  at = 1, horizontal = T, add = TRUE, xaxt = "n", yaxt = "n", xlab = '', ylab = '', width = 10, pch = 19, cex = .3, boxwex = 4, lwd = .4, col = "gold", border = "grey");
#  text(mean(dta_box$trait), 5, 'Days to anthesis', cex = .5)
 }
 dev.off()

}
#####
Plot_Trait_mean_envParas <- function( env_mean_trait, env_paras, d1, d2, trait, exp_dir, env_cols){
  days <- c(d1:d2); 
  env_facts_matrix <- matrix(nrow = nrow(env_mean_trait), ncol = 6);
  for (e_i in 1:nrow(env_mean_trait)) {
    e <- env_mean_trait$env_code[e_i];
    env_para <- subset(env_paras, env_paras$env_code == e);
    env_mean <- colMeans(env_para[days, c(8, 9, 10, 11, 12) - 3]); ### DL, GDD, PTT, PTR, PTS
    
    env_facts_matrix[e_i,] <- c(env_mean_trait$meanY[e_i], round(env_mean, 4) );
  }
  colnames(env_facts_matrix) <- c('env_code', 'DL', 'GDD', 'PTT', 'PTR', 'PTS');
#  trait_mean_envPara <- merge(env_mean_trait, env_facts_matrix, stringsAsFactors = F);
  
  pdf_file <- paste(exp_trait_dir, trait, 'Mean_', nrow(env_mean_trait), 'EnvPara.pdf', sep = ''); 
  pdf(pdf_file,width= 7.5,height= 1.5,pointsize=6)

  layout(matrix(c(1:5), ncol = 5));
  
  for (i in 1:5) {
   par(mar = c(2.5, 2.0, 1, 0.5) , mgp = c(0.7, 0.01, 0), tck = -0.01, family = "mono");
   plot(env_facts_matrix[, i + 1], env_facts_matrix[,1], xlab = colnames(env_paras)[i + 4], ylab = paste(trait, ' mean', sep = ''),  pch = 19, col = env_cols);
   abline(lm(env_facts_matrix[,1] ~ env_facts_matrix[, i + 1]), lty = 2);
   legend("bottom", paste('r = ', round(cor(env_facts_matrix[,1] , env_facts_matrix[, i + 1]), 3), sep = ''), bty = "n")
   if (i == 1) { legend("topleft", env_mean_trait$env_code, pch = 19, col = env_cols, bty = "n" )};
  }
  dev.off()
  
 
 

}

#### convert FT in julian day to days-after-planting and GDD
Convert_HDj_HDdap_HDgdd <- function(exp_traits, PTT_PTR, exp_traits_file, FTj_tag) {
  env_codes <- unique(exp_traits$env_code);
  for (e_i in 1:length(env_codes)) {
    e <- env_codes[e_i];
    
    env_ptt_ptr <- subset(PTT_PTR, PTT_PTR$env_code == e);
    gdds <- cumsum(env_ptt_ptr$GDD);
    
    if (FTj_tag == 1) {
     env_trait <- subset(exp_traits, exp_traits$env_code == e & !is.na(exp_traits$FTj));
     dop <- env_ptt_ptr$date[1];
     doy <- as.numeric(strftime(dop, format = "%j"));
     ftdap <- floor(env_trait$FTj - doy + 1);
     ftdap[ftdap <0 ] <- ftdap[ftdap < 0] + 365
    } else {
     env_trait <- subset(exp_traits, exp_traits$env_code == e & !is.na(exp_traits$FTdap));
     ftdap <- env_trait$FTdap;
    }
    ftgdd <- gdds[ftdap];
    t_df <- data.frame(line_code = env_trait$line_code, env_code = env_trait$env_code, FTdap = ftdap, FTgdd = ftgdd);
    if (e_i == 1) { FT_df <- t_df} else { FT_df <- rbind(FT_df, t_df)};
  }
  exp_traits <- merge(exp_traits, FT_df, all.x = T);
  write.table(exp_traits, file = paste(exp_traits_file, '_addFTgdd', sep = '' ), sep = "\t", quote = F, row.name = F)
}



### #### stations from ISD
#local_isd_st_file <- paste(Dir, 'all_isd_stations', sep = '');
#local_isd_target_st_file <- paste(sp_isd_dir, '0target_isd_stations', sep = '');
#if (!file.exists(local_isd_target_st_file)) {
#  if (!file.exists(local_isd_st_file)) { isd_all_sts <- isd_stations(); write.csv(isd_all_sts, local_isd_st_file)} else {isd_all_sts <- read.csv(local_isd_st_file);}
#  sts_target <- dplyr::filter(isd_all_sts, begin <= exp_s_day & end >= exp_e_day 
#                                          & between(lat,lat_range[1] - 2, lat_range[2] + 2) 
#                                          & between(lon,lon_range[1] - 2, lon_range[2] + 2)
#                                          );
#  write.csv(sts_target, local_isd_target_st_file  )                                          
#}
#STS_isd <- read.csv(local_isd_target_st_file);      
