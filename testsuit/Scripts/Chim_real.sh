#!/bin/bash
grep 'Sum Chi_m :' out|awk '{print $4}'
