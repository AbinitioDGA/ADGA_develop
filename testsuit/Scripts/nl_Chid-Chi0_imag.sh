#!/bin/bash
grep 'Sum Chi_d^nl - Chi_0^q:' output-test/out|awk '{print $6}'
