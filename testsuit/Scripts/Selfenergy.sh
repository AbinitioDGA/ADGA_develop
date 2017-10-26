#!/bin/bash
grep 'Tr[Self-energy] :' out|awk '{print $3}'
