#!/bin/bash
grep 'Sum Chi_d^q:' output-test/out|awk '{print $3}'
