#!/bin/bash

cat $1 | grep -v "^#" | cut -f 3 | sort | uniq -c
