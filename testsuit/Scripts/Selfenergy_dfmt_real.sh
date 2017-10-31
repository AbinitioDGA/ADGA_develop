#!/bin/bash
grep 'Tr\[Local Self-energy\] :' out|awk '{print $3}'
