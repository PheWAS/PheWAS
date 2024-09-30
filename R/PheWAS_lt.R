# 
# phewas_lt <- function(demos, pheno, phecode_info = PheWASmaps::pheinfo){
# # #setnames(demos, 'GENDER', 'gender')
# # head(pheno)
# # pheno <- as.data.table(pheno)
#  pheno=pheno[GRID %in% demos$GRID] #exclude GRIDs in pheno that are not in demos
# # #setnames(pheno, "n", "N")
# # 
# # #phecode_info <- fread("c5_x_prev 9-13-2023.csv") #I verified this file's phecode_sex matches BP_PheWAS_exp_7_8_phecodeX 1-22-2024
# # is.data.table(phecode_info)
# # phecode_info <- as.data.table(phecode_info)
# # head(phecode_info)
# # #phecode_info <- phecode_info[, .(phecode, phecode_string, phecode_sex = sex)]
# # head(phecode_info)
# # #manual changes to incorporate changes Lisa has made to phecodeX
# # phecode_info[phecode == "GU_603.6", phecode_sex := "Male"] #Scrotal pain*
# # phecode_info[phecode == "GU_612.2", phecode_sex := "Both"] #Breast conditions, congenital or relating to hormones
# # phecode_info[phecode == "SS_834,1", phecode_sex := "Female"] #Abnormal Papanicolaou smear
# # 
# # #var is the categorical variable, code requires a factor with values 10, 11, 12
# # setnames(demos, "ps_n", "var")
# # 
# # # set case selection criteria, phecode_sex_method, & initialize result data tables 
# # #phecode_sex method 1 = empirical, 2 = data derived
#  ## set as input phecode_sex_method <- "1"
# # 
# # #covariates 1=age & sex, 2 = age & sex & PP
# # covariates <- "1"
# # 
# # #predictor p11, p12, or both
# # predictor <- "both"
# # 
# # min_cases_to_run = 400
# # cur_MCC = 2
# # 
# # results=data.table()
# # result=data.table()
# # 
# # # make a summary of phecodes and the number of available cases (MCC >= cur_MCC)
# # pheno_summary = pheno[N >= cur_MCC, .N, by='phecode']
# # phecodes_to_run = pheno_summary[N >= min_cases_to_run]$phecode
# # length(phecodes_to_run) #c5 = 1857
# # j = 3
# # 
# # 
# for(j in 1:length(phecodes_to_run)){
# 
#   cur_phecode = phecodes_to_run[j]
#   phecode_string = phecode_info[phecode == cur_phecode]$phecode_string
#   print('[' %++% j %++% ']  ' %++% Sys.time() %++%  '--' %++% cur_phecode %++% ' -- ' %++% phecode_string)
# 
#   d=get_pheno(pheno,demos,phecode=cur_phecode, MCC=cur_MCC)
# 
#   # check phecode sex
#   if(phecode_sex_method == 1){
#     head(d)
#     # empirical phecode_sex
#     cur_sex = phecode_info[phecode == cur_phecode]$phecode_sex
#     if(cur_sex == 'Male' ){
#       d = d[gender == 'M']
#     } else if(cur_sex == 'Female' ){
#       d = d[gender == 'F']
#     }
#   } else {
# 
#     # data-derived phecode_sex
#     cur_sex = 'Both'
#     if(nrow(d[gender == 'M' & pheno == 1]) < nrow(d[pheno == 1])/10){
#       cur_sex <- 'Female'
#       d = d[gender == 'F']
#     } else if(nrow(d[gender == 'F' & pheno ==1]) < nrow(d[pheno == 1])/10){
#       cur_sex <- 'Male'
#       d <- d[gender == 'M']
#     }
#   }
# 
#   cases=table(d$pheno)[["1"]]
#   controls=table(d$pheno)[["0"]]
#   missing = nrow(d[is.na(pheno)])
# 
#   ## exclude GRIDs with missing pheno
#   d = d[!is.na(pheno)]
# 
#   p10 = nrow(d[var == '10'])
#   p10_case = nrow(d[var=='10' & pheno == 1])
#   p10_control = nrow(d[var=='10' & pheno == 0])
# 
#   if(cur_sex == 'Both'){
#     if(covariates == 1){
#       m=glm(pheno~var + last_age + gender, data=d, family="binomial")
#     } else {
#       m=glm(pheno~var + mdn_pp + last_age + gender, data=d, family="binomial")
#     }
#   } else {
#     if(covariates == 1){
#       m=glm(pheno~var + last_age, data=d, family="binomial")
#     } else {
#       m=glm(pheno~var + mdn_pp + last_age, data=d, family="binomial")
#     }
#   }
# 
#   #summary(m)
#   #plot(allEffects(m))
# 
#   #extract model parameters
#   conf=confint.default(m)
#   coefficients(summary(m))
# 
#   #use car::vif: to generate model variance inflation factor
#   vif <- vif(m)
#   print(m)
#   vif_var <- vif['var', 'GVIF'] #extract GVIF for 'var'
#   print(vif)
#   if(predictor == 'p11' | predictor == 'both'){
# 
#     p11 = nrow(d[var == '11'])
#     p11_case = nrow(d[var=='11' & pheno == 1])
#     p11_control = nrow(d[var=='11' & pheno == 0])
#     P_p11 = coefficients(summary(m))["varp11","Pr(>|z|)"]
#     OR_p11 = exp(coefficients(summary(m))["varp11","Estimate"])
#     L95_p11 = exp(conf['varp11',1])
#     U95_p11 = exp(conf['varp11',2])
# 
#     result_p11 = data.table(
#       phecode = cur_phecode, phecode_string = phecode_string, phecode_sex = cur_sex, BP_category = 'p11',
#       cases = cases, controls = controls, n_exclude = missing,
#       P = P_p11, OR = OR_p11, L95 = L95_p11, U95 = U95_p11,
#       aff = p11, unaff = p10,
#       aff_case = p11_case, aff_control = p11_control,
#       unaff_case = p10_case, unaff_control = p10_control,
#       vif_var = vif_var
#     )
# 
#     results=rbind(results, result_p11)
#   }
# 
#   if(predictor == 'p12' | predictor == 'both'){
# 
#     p12 = nrow(d[var == '12'])
#     p12_case = nrow(d[var=='12' & pheno == 1])
#     p12_control = nrow(d[var=='12' & pheno == 0])
# 
#     P_p12 = coefficients(summary(m))["varp12","Pr(>|z|)"]
#     OR_p12 =exp (coefficients(summary(m))["varp12","Estimate"])
#     L95_p12 = exp(conf['varp12',1])
#     U95_p12 = exp(conf['varp12',2])
# 
#     result_p12 = data.table(
#       phecode = cur_phecode, phecode_string = phecode_string, phecode_sex = cur_sex, BP_category = 'p12',
#       cases = cases, controls = controls, n_exclude = missing,
#       P = P_p12, OR = OR_p12, L95 = L95_p12, U95 = U95_p12,
#       aff = p12, unaff = p10,
#       aff_case = p12_case, aff_control = p12_control,
#       unaff_case = p10_case, unaff_control = p10_control,
#       vif_var = vif_var
#     )
# 
#     results=rbind(results, result_p12)
#   }
# }
# }
