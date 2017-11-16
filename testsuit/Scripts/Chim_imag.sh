#!/bin/bash
grep 'Sum Chi_m^q:' output-test/out|awk '{print $4}'
