#!/bin/bash
grep 'Sum Chi_m^nl - Chi_0^q:' output-test/out|awk '{print $5}'
