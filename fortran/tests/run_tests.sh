#!/bin/bash
for LANGUAGE in F77 F95; do
  echo "Running ${LANGUAGE} tests"
  make -s test LANGUAGE=${LANGUAGE}
done
