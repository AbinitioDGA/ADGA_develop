#!/bin/bash
grep 'Sum Chi_m^nl:' output-test/out|awk '{print $4}'
