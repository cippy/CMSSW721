#!/bin/bash
echo
./tmp/main config/monojet_signalRegion_config_spring15_25ns.txt
echo
echo
./tmp/main config/zmumujets_ControlRegion_config_spring15_25ns.txt
echo
echo
./tmp/main config/wmunujets_ControlRegion_config_spring15_25ns.txt
echo
echo
#./tmp/main config/zeejets_ControlRegion_config_spring15_25ns.txt
#echo
#echo
#./tmp/main config/zeejets_ControlRegion_config_spring15_25ns.txt -calibEle
#echo
#echo
#./tmp/main config/wenujets_ControlRegion_config_spring15_25ns.txt
#echo
#echo
#./tmp/main config/wenujets_ControlRegion_config_spring15_25ns.txt -calibEle
#echo
#echo
echo "producing plots for monojet (exclusive)"
./macro/distribution 0
echo
echo
echo "producing plots for monoV"
./macro/distribution 1
echo
echo
echo "========================"
echo "===  End of Analysis ==="
echo "========================"
echo "Checking for running jobs"
bjobs
echo
