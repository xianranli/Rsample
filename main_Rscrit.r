
library(rnoaa);
library(dplyr);
library(RCurl);
library(htmltab);
library(corrgram);
library(lubridate);
library(animation);
library("colorspace");
library(RColorBrewer);
##
source('D:/0GbE/scripts/sub_funcs.r');
######################################
col_wdw <- 25;
col_palette <- diverge_hcl(col_wdw + 1, h = c(260, 0), c = 100, l = c(50, 90), power = 1)

Dir <- 'D:/0GbE/';
sp <- 'Oat';
exps <- 'CoreS';

exp_s_year <- 2010; #exp_s_day <- paste(exp_s_year, '0101', sep = "\t"); 
exp_e_year <- 2011; #exp_e_day <- paste(exp_e_year, '1231', sep = "\t");

sp_dir <- paste(Dir, sp, '/', sep = '');
exp_dir <- paste(sp_dir, exps, '/', sep = '');

## for barley and wheat
t_base <- 32; t_max1 <- 95; t_max2 <- 7000; Haun_threshold <- 0;

## for barley and wheat
#t_base <- 32; t_max1 <- 95; t_max2 <- 70; Haun_threshold <- 384;
## for maize and sorghum
# t_base <- 50; t_max1 <- 100; t_max2 <- 1000; Haun_threshold <- 0;

env_meta_file <- paste(exp_dir, '12Env_meta_table', sep = ''); ## S2_Met_22Env_meta_table
env_meta_info_0 <- read.table(env_meta_file, header = T, sep = "\t", stringsAsFactors = F);

searching_daps <- 150;
ptt_ptr_file <- paste(exp_dir, nrow(env_meta_info_0), 'Envs_PTT_PTR_DAP',  searching_daps, sep = ''); ##  '_adjGDD',
if (!file.exists(ptt_ptr_file)) { Compile_PTT_PTR(exp_dir, env_meta_info_0, exp_s_year, exp_e_year, searching_daps, ptt_ptr_file) };
PTT_PTR <- read.table(ptt_ptr_file, header = T , sep = "\t");

exp_traits_file <- paste(exp_dir, 'traits_ori_addFT', sep = '');
exp_traits <- read.table(exp_traits_file, sep = "\t", header = T, stringsAsFactors = F);
#exp_traits <- subset(exp_traits, exp_traits$gen == 0);

#FTj_tag <- 0;  Convert_HDj_HDdap_HDgdd(exp_traits, PTT_PTR, exp_traits_file, FTj_tag); #

#exp_traits_file <- paste(exp_dir, 'paired_traits_ori_trim_addFT', sep = '');
#exp_traits <- read.table(exp_traits_file, sep = "\t", header = T, stringsAsFactors = F);
#exp_traits <- exp_traits[-grep('NIL', exp_traits$line_code),]
all_env_codes <- unique(exp_traits$env_code);
env_cols <- rainbow_hcl(length(all_env_codes), c = 80, l = 60, start = 0, end = 300, fixup = TRUE, alpha = 0.75);

trait <- 'GY'; 
exp_trait_dir <- paste(exp_dir, trait, '/', sep = ''); if (!dir.exists(exp_trait_dir))  { dir.create(exp_trait_dir)};

lInd <- which(colnames(exp_traits) == 'line_code'); eInd <- which(colnames(exp_traits) == 'env_code'); tInd <- which(colnames(exp_traits) == trait);

exp_trait <- exp_traits[,c(lInd, eInd, tInd)]; ### make sure the colname is line env trait
#line_outliers <- c('Class1', 'Class2', 'Class3a', 'Class3b', 'Class3c', 'Class4');
#exp_trait <- exp_trait[!(exp_traits$line_code %in% line_outliers),]
colnames(exp_trait)[3] <- 'Yobs';
exp_trait <- exp_trait[!is.na(exp_trait$Yobs),];

#aggregate(Yobs ~ env_code, data = exp_trait, function(x) sum( !is.na(x) ))
#aggregate(Yobs ~ line_code, data = exp_trait, function(x) sum( !is.na(x) ))

### pairwise correlations among trials; trait distribution among trials
Pairwise_trait_env_distribution_plot(exp_trait, exp_trait_dir, trait, all_env_codes, env_meta_info_0);
## PH
#env_outliers <- c('U06_WestLafayette'); ##c('UEOPN_2000_Brookings', 'UMOPN_2007_Morris', 'UMOPN_2005_Watertown', 'UMOPN_2007_Watertown'); ##, 'HNY15', 'KNY15', 'BZD15', 'CRM15', , 'CHM16'
## FTgdd
#env_outliers <- c('CS10_Tetonia');
##line_outliers <- c('WINONA');
#exp_trait <- exp_trait[!(exp_trait$env_code %in% env_outliers),];
#exp_trait <- exp_trait[!(exp_trait$line_code %in% line_outliers),];

line_codes <- unique(exp_trait$line_code); 
env_mean_trait_0 <- na.omit(aggregate(x = exp_trait$Yobs, by = list(env_code = exp_trait$env_code), mean, na.rm = T));
colnames(env_mean_trait_0)[2] <- 'meanY';
env_mean_trait <- env_mean_trait_0[order(env_mean_trait_0$meanY),];
##### searching the critical window with the highest correlation with environment mean
##### the window can be adjusted based on biological meaning
#searching_daps <- 150;
Exhaustive_search(env_mean_trait, PTT_PTR, searching_daps, exp_trait_dir, exp_traits$FTdap, trait, 1, searching_daps, searching_daps)#; searching_daps, searching_daps);

#### plot trait mean and environment paramenters
maxR_dap1 <- 51;
maxR_dap2 <- 78;
Plot_Trait_mean_envParas(env_mean_trait, PTT_PTR, maxR_dap1, maxR_dap2, trait, exp_trait_dir, env_cols);

#### Leave-one-environment-out cross validation
p <- 1;
maxR_dap1 <- '96'; maxR_dap2 <- '114'; PTT_PTR_ind <- 8; ## PTT -> 7; PTR -> 8;
obs_prd_file <- paste(exp_trait_dir, trait, '_', nrow(env_mean_trait), 'Env_LOO_by_Lines_kPara_D', maxR_dap1, '_', maxR_dap2, sep = '');
LOO_pdf_file <- paste(exp_trait_dir, trait, '_', nrow(env_mean_trait), 'Env_LOO_by_Lines_kPara_D', maxR_dap1, '_', maxR_dap2, '.pdf', sep = '');
if (!file.exists(obs_prd_file)) {
 
 maxR_win <- c(maxR_dap1:maxR_dap2);
 prdM <- env_mean_trait;
 maxR_envPara <- matrix(ncol = 2, nrow = nrow(env_mean_trait));
 kPara <- c();
 for (e_i in 1:nrow(env_mean_trait)) {
   envParas <- subset(PTT_PTR, PTT_PTR$env_code == env_mean_trait$env_code[e_i]);
   envPara <- mean(envParas[maxR_win, PTT_PTR_ind]);
#   envPara <- mean(envParas$PTR[maxR_win]);
   kPara[e_i] <- envPara;
 }
 prdM$kPara <- kPara;
 obs_prd_m <- matrix(0, ncol = 7,nrow = nrow(exp_trait));
 
 n <- 0; 
 
 for (l in line_codes) {
   l_trait <- subset(exp_trait, exp_trait$line == l); l_mean <- mean(l_trait$Yobs, na.rm = T);
   ril_data <- merge(prdM, l_trait,  all.x = T);
     if (length(which(!is.na(ril_data$Yobs))) > 3) {
     for (e in 1:nrow(ril_data)) {
       obs_trait <- ril_data$Yobs[e];
       if (!is.na(obs_trait)) {
         trn <- ril_data[-e,];
         prd_trait_mean  <- round(predict( lm(Yobs ~ meanY, data = trn), ril_data[e,]), digit = 3);
         prd_trait_kpara <- round(predict( lm(Yobs ~ kPara, data = trn), ril_data[e,]), digit = 3);
         n <- n + 1;
         obs_prd_m[n,] <- c(ril_data$env_code[e], p, l, prd_trait_mean, prd_trait_kpara, obs_trait, l_mean);
         
       }
     }
   }
 }
 
 obs_prd_m <- obs_prd_m[1:n,]
 colnames(obs_prd_m) <- c('env_code', 'pop_code', 'ril_code', 'Prd_trait_mean', 'Prd_trait_kPara', 'Obs_trait', 'Line_mean');
 write.table(obs_prd_m, file = obs_prd_file, sep = "\t", quote = F);

}

Obs_Prd_m <- read.table(obs_prd_file, sep = "\t", header = T);
env_rs <- matrix(ncol = 3, nrow = nrow(env_mean_trait));
for (e_i in 1:nrow(env_mean_trait)) {
  env_obs_prd <- subset(Obs_Prd_m, Obs_Prd_m$env_code == env_mean_trait$env_code[e_i]);
  env_rs[e_i,] <- c( sprintf( "%.2f", cor(env_obs_prd[,4], env_obs_prd[,6])), sprintf( "%.2f", cor(env_obs_prd[,5], env_obs_prd[,6])), sprintf( "%.2f", cor(env_obs_prd[,7], env_obs_prd[,6]))); 
}

xy_lim <- range(Obs_Prd_m[4:6])

pdf(LOO_pdf_file ,width = 4, height= 4,pointsize=6)
#
#for (p in Pops) {
layout(matrix(c(1:4), 2, 2, byrow = T));
# screen(1)
 obs_prd_m <- subset(Obs_Prd_m, Obs_Prd_m$pop_code == p);
 par(mar = c(2.5, 2.5, 1.0, 0.5) , mgp = c(1, 0.25, 0), tck = -0.005, family = "mono");
 plot(obs_prd_m[,6], obs_prd_m[,4], col = env_cols[match(obs_prd_m[,1], all_env_codes)], pch = 19, cex = .4, ylab = paste('Predicted ', trait, ' by envMean', sep = ''), xlab = paste('Observed ', trait, '', sep = ''), xlim = xy_lim, ylim = xy_lim);
 abline(a = 0, b = 1, lty = 2, col = "gray59");
 r1 <- sprintf("%.2f", cor(obs_prd_m[,6], obs_prd_m[,4]));
 legend("top", legend= substitute(paste(italic('r'), " = ", R1), list(R1 = r1)), bty = "n")
 LGs <- c()
 for (e in 1:nrow(env_mean_trait)) {
   LGs <- append(LGs, bquote(.(E) * ' (' *italic(r) == .(A) * ')', list(E = as.vector(env_mean_trait$env_code)[e], A = env_rs[e,1])))
 }
 legend("topleft", legend=do.call("expression", LGs), col = env_cols[match(env_mean_trait$env_code, all_env_codes)],  pch = 19, cex = .65, bty = "n")

 mtext('A', side = 3, at = xy_lim[1], line = 0.1, cex = .8);
# screen(2)
 par(mar = c(2.5, 2.5, 1.0, 0.5) , mgp = c(1, 0.25, 0), tck = -0.005, family = "mono");
 plot(obs_prd_m[,6], obs_prd_m[,5], col = env_cols[match(obs_prd_m[,1], all_env_codes)], pch = 19, cex = .4, ylab = paste('Predicted ', trait, ' by envPara', sep = ''), xlab = paste('Observed ', trait, '', sep = ''), xlim = xy_lim, ylim = xy_lim);
 abline(a = 0, b = 1, lty = 2, col = "gray59");
 mtext('B', side = 3, at = xy_lim[1], line = 0.1, cex = .8);
 r2 <- sprintf("%.2f", cor(obs_prd_m[,6], obs_prd_m[,5]));
 legend("top", legend= substitute(paste(italic('r'), " = ", R1), list(R1 = r2)), bty = "n")
 LGs <- c()
 for (e in 1:nrow(env_mean_trait)) {
   LGs <- append(LGs, bquote(.(E) * ' (' *italic(r) == .(A) * ')', list(E = as.vector(env_mean_trait$env_code)[e], A = env_rs[e,2])))
 }
 legend("topleft", legend=do.call("expression", LGs), col = env_cols[match(env_mean_trait$env_code, all_env_codes)],  pch = 19, cex = .65, bty = "n")

#}
 par(mar = c(2.5, 2.5, 1.0, 0.5) , mgp = c(1, 0.25, 0), tck = -0.005, family = "mono");
 plot(prdM$kPara, prdM$meanY, col = env_cols[match(prdM$env_code, all_env_codes)],  ylim = xy_lim, pch = 19, cex = .65, xlab = 'envPara', ylab = 'Observed population mean');
 mtext(prdM$env_code, side = 1, at = prdM$kPara, las = 2, line = -2, cex = .6 )
 abline(lm(prdM$meanY ~ prdM$kPara))
  r2 <- sprintf("%.2f", cor(prdM$meanY, prdM$kPara));
 legend("top", legend= substitute(paste(italic('r'), " = ", R1), list(R1 = r2)), bty = "n")

 mtext('C', side = 3, at = min(prdM$kPara), line = 0.1, cex = .8);


# screen(3)
 par(mar = c(2.5, 2.5, 1.0, 0.5) , mgp = c(1, 0.25, 0), tck = -0.005, family = "mono");
 plot(obs_prd_m[,6], obs_prd_m[,7], col = env_cols[match(obs_prd_m[,1], all_env_codes)], pch = 19, cex = .4, xlab = paste('Observed ', trait, sep = ''), ylab = paste('Predicted ', trait, ' by BLUE', sep = ''), xlim = xy_lim, ylim = xy_lim);
 abline(a = 0, b = 1, lty = 2, col = "gray59");
 r1 <- sprintf("%.2f", cor(obs_prd_m[,6], obs_prd_m[,7]));
 legend("top", legend= substitute(paste(italic('r'), " = ", R1), list(R1 = r1)), bty = "n")
 LGs <- c()
 for (e in 1:nrow(env_mean_trait)) {
   LGs <- append(LGs, bquote(.(E) * ' (' *italic(r) == .(A) * ')', list(E = as.vector(env_mean_trait$env_code)[e], A = env_rs[e,3])))
 }
 legend("topleft", legend=do.call("expression", LGs), col = env_cols[match(env_mean_trait$env_code, all_env_codes)],  pch = 19, cex = .65, bty = "n")
 mtext('D', side = 3, at = xy_lim[1], line = 0.1, cex = .8);


dev.off()

