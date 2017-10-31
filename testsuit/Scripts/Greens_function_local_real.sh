#!/bin/bash
grep 'Tr\[Local Greens function\] :' out|awk '{print $3}'
