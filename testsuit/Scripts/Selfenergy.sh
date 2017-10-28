#!/bin/bash
grep 'Tr\[Self-energy\] :' output_test/out|awk '{print $3}'
