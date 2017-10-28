#!/bin/bash
grep 'Tr\[Chi_m\] :' out|awk '{print $3}'
