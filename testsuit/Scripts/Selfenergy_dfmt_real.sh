#!/bin/bash
grep 'Tr\[Local Self-energy\]:' output-test/out|awk '{print $3}'
