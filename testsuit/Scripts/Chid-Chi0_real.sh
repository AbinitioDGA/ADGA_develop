#!/bin/bash
grep 'Sum Chi_d^q  - Chi_0^q:' output-test/out|awk '{print $5}'
