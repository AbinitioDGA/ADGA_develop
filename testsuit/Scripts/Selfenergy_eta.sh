#!/bin/bash
grep 'Tr\[Self-energy\[eta\]\] :' out|awk '{print $3}'
