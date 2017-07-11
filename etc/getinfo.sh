#!/bin/bash

#hardware information
lshw -short > info.txt

#detailed cpu info
lscpu >> info.txt

#distribution and kernel version
uname -a >> info.txt
