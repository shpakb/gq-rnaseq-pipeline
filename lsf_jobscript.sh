#!/bin/bash
# invoke with --jobscript

# Propagate LSF tmp dir:
# TMPDIR for LSF nodes
if [ -n "$__LSF_JOB_TMPDIR__" ]; then
  export TMPDIR=$__LSF_JOB_TMPDIR__
fi

# properties = {properties}
{exec_job}