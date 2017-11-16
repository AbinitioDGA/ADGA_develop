#!/bin/bash
grep 'Sum Chi_0^w:' output-test/out|awk '{print $3}'
