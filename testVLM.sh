#!/bin/sh

julia unit_test_run_VLM.jl
diff check_output.csv original_output.csv
