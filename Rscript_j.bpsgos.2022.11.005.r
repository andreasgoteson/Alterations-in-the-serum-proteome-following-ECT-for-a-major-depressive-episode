### Code for gls-analyses from Göteson et al, BP:GOS 2022
### https://doi.org/10.1016/j.bpsgos.2022.11.005

### NOTE! This is example code and will not run without access to original files
### The file is a long-format data frame nested by protein

# Define models
run_gls = function(df, model="primary"){
  
  if (model=="primary")
    model = with(df, gls(NPX ~ t + AGE_AT_SAMPLE + SEX + PlateID, 
                       corr=corSymm(form= ~ time | STUDYPERSONID), 
                       weights = varIdent(form = ~ 1 | t)))
  
  else if (model=="secondary")
    model = with(df, gls(NPX ~ t + AGE_AT_SAMPLE + SEX + PlateID + LITIUM_BEFORE + neuroleptics + tobac1 + Q20_LADDNINGFORSTA, 
                       corr=corSymm(form= ~ time | STUDYPERSONID), 
                       weights = varIdent(form = ~ 1 | t)))
  
  else if (model=="interaction")
    model = with(df, gls(NPX ~ response*t + AGE_AT_SAMPLE + SEX + PlateID, 
                       corr=corSymm(form= ~ time | STUDYPERSONID), 
                       weights = varIdent(form = ~ 1 | t)))
    
  else
    print("no valid model selected")

  return(model)
}

# And a wrapper to run the function on nested data and extract statistics
gls_wrapper = function(nested_data, model="primary"){
  
  res = list()
  for (i in 1:nrow(nested_data)){
    res[[i]] = run_gls(nested_data[i,]$data[[1]], 
                       model = model)
  }
  
  names(res) = nested_data$Assay
  
  bind_res = lapply(res, function(x){
    x %>% summary %>% coefficients %>% 
      data.frame %>% rownames_to_column(var="term")
    }) %>% bind_rows(.id = "Assay") %>% 
    # Add FDR
    group_by(term) %>% 
    mutate(FDR = p.adjust(p.value, method="BH")) %>% 
    ungroup
  
  # Return a list of the models, data, and bind_res
  out = list(data = nested_data,
             models = res,
             results = bind_res)
}

# Now create an object with nested data and run the analyses:
nested_data = prefect %>%  
  dplyr::select(Assay, NPX, t, STUDYPERSONID, 
                AGE_AT_SAMPLE, SEX, PlateID) %>% 
  # Prepare data
  mutate(time=as.integer(t)) %>% 
  arrange(t) %>% 
  group_by(Assay) %>% nest(data = -Assay)

# Run models
out1 = gls_wrapper(nested_data, model="primary")