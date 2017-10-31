#!/bin/bash
grep 'Tr\[Total local Self-energy\]:' output-test/out|awk '{print $5}'
