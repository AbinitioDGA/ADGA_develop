#!/bin/bash
grep 'Tr\[Total loc Self-energy\] :' output_test/out|awk '{print $4}'
