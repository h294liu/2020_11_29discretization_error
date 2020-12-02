#Configuration File for Ostrich Program
ProgramType PADDS
OstrichWarmStart no

BeginFilePairs
step2_generate_HRU_tpl.py; model/step2_generate_HRU.py
EndFilePairs

ObjectiveFunction 	GCOP
ModelExecutable     run_model.sh

#Parameter Specification
BeginParams
  #parameter	init.	low	high tx_in tx_ost tx_out  format
	slp_thrsh1_value	random	0.001	52.0794	none	none	none free
	slp_thrsh2_value	random	0.001	52.0794	none	none	none free
	slp_thrsh3_value	random	0.001	52.0794	none	none	none free
	asp_thrsh1_value	random	0.001	179.974	none	none	none free
	asp_thrsh2_value	random	0.001	179.974	none	none	none free
	asp_thrsh3_value	random	0.001	179.974	none	none	none free
EndParams

BeginResponseVars
  #name      filename                          keyword     line col  token
  HRU_num      ROOT_DIR/CASE/Diagnostics.txt;  OST_NULL   1    1    ','
  Sw_Rad_error        ROOT_DIR/CASE/Diagnostics.txt;  OST_NULL   1    2    ','
EndResponseVars

BeginGCOP
  CostFunction HRU_num 
  CostFunction Sw_Rad_error 
  PenaltyFunction APM
EndGCOP

BeginConstraints
  # name     type     penalty    lwr   upr   resp.var
EndConstraints

BeginPADDS
  MaxIterations MAXITERATIONS
  PerturbationValue 0.2
#  SelectionMetric Random
#  SelectionMetric CrowdingDistance
#  SelectionMetric EstimatedHyperVolumeContribution
  SelectionMetric ExactHyperVolumeContribution
EndPADDS