#!/bin/bash
grep 'Tr\[Local Greens function\]:' output-test/out|awk '{print $5}'
