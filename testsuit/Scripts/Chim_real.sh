#!/bin/bash
grep 'Sum Chi_m:' output-test/out|awk '{print $3}'
