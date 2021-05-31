#!/bin/sh
/Users/ekellogg/src/Rosetta/main/source/bin/rosetta_scripts.default.macosclangrelease \
-parser:protocol $1 \
-database /Users/ekellogg/src/Rosetta/main/database \
-edensity:mapfile $3 \
-s $2 \
-default_max_cycles 200 \
-ignore_unrecognized_res \
-nstruct 1 \
-overwrite \
-parser::script_vars readbeams=$4 beams=$5 steps=$6 pcount=$7 \
-parser::script_vars filterprevious=$8 filterbeams=$9 \
-edensity:sliding_window 3 \
-mapreso 3.2 \
-corrections::beta_nov16 \
-missing_density_to_jump
