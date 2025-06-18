#!/bin/bash

#
# Script used to generate Extended CFF files (test_data/cff_extended) from testing datasets (test_data/caller_output_files)
#

# BEERS
/MetaFusion/cff_generation/RUN_convert_fusion_results_to_cff.sh -c /MetaFusion/test_data/caller_output_files -s /MetaFusion/test_data/caller_output_files/sample_info.txt -d BEERS -o /MetaFusion/test_data/cff_extended/BEERS -t "arriba defuse ericscript fusionmap integrate star_fusion star_seqr"

# BRCA
/MetaFusion/cff_generation/RUN_convert_fusion_results_to_cff.sh -c /MetaFusion/test_data/caller_output_files -s /MetaFusion/test_data/caller_output_files/sample_info.txt -d BRCA -o /MetaFusion/test_data/cff_extended/BRCA -t "arriba defuse ericscript fusionmap integrate star_fusion star_seqr"

# DREAM
/MetaFusion/cff_generation/RUN_convert_fusion_results_to_cff.sh -c /MetaFusion/test_data/caller_output_files -s /MetaFusion/test_data/caller_output_files/sample_info.txt -d DREAM.SIM45.SIM52 -o /MetaFusion/test_data/cff_extended/DREAM.SIM45.SIM52 -t "arriba defuse ericscript fusionmap integrate star_fusion star_seqr"

# MELANOMA
/MetaFusion/cff_generation/RUN_convert_fusion_results_to_cff.sh -c /MetaFusion/test_data/caller_output_files -s /MetaFusion/test_data/caller_output_files/sample_info.txt -d MELANOMA -o /MetaFusion/test_data/cff_extended/MELANOMA -t "arriba defuse ericscript fusionmap integrate star_fusion star_seqr"

# SIM50
/MetaFusion/cff_generation/RUN_convert_fusion_results_to_cff.sh -c /MetaFusion/test_data/caller_output_files -s /MetaFusion/test_data/caller_output_files/sample_info.txt -d SIM50 -o /MetaFusion/test_data/cff_extended/SIM50 -t "arriba defuse ericscript fusionmap integrate star_fusion star_seqr"

# SIM101
/MetaFusion/cff_generation/RUN_convert_fusion_results_to_cff.sh -c /MetaFusion/test_data/caller_output_files -s /MetaFusion/test_data/caller_output_files/sample_info.txt -d SIM101 -o /MetaFusion/test_data/cff_extended/SIM101 -t "arriba defuse ericscript fusionmap integrate star_fusion star_seqr"