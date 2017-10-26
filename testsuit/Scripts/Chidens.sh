#!/bin/bash
grep 'Tr[Chi_d] :' out|awk '{print $3}'
