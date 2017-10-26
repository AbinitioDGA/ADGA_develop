#!/bin/bash
grep 'Tr[Self-energy[gamma]] :' out|awk '{print $3}'
