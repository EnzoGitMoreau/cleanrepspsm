#!/bin/zsh

# Loop through numbers 1 to 8
for i in {1..13}
do
  # Run the program with the current number as an argument
  ./tests $(5*$i) 12000 100 0.5

done
