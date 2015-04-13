addPhewasGroups = function(d, keep_newlines=T) {
  if(!length(d$phenotype)) stop("d must contain attribute 'phenotype'.")
  if(!keep_newlines) code_group_ranges$group=sub("\\n"," ",code_group_ranges$group)
  d=merge(d,code_group_ranges,by=c())
  d=d[d$phenotype<d$max & d$phenotype>=d$min ,!(names(d) %in% c("min","max"))]
  #If there are any groups with no phenotypes, warn the user
  present_groups=code_group_ranges$groupnum %in% unique(d$groupnum)
  missing_groups=code_group_ranges$range_str[!(present_groups)]
  if(length(missing_groups)>0) {
    warning("No codes found for group(s) " %+% paste0(missing_groups,collapse=", ") %+% ".")
  }
  d
}