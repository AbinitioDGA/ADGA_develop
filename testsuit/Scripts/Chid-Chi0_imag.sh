#!/bin/bash
grep 'Sum Chi_d - Chi_0^q :' out|awk '{print $7}'
