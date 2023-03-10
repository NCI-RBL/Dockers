#!/usr/bin/env bash
### Script to get N random lines from a file
# Requires coreutils (brew install coreutils)

input_file=  number=

while getopts i:n: opt; do
  case $opt in
  i)
      input_file=$OPTARG
      ;;
  n)
      number=$OPTARG
      ;;
  \?)
      echo "Invalid option: -$OPTARG" >&2
      ;;
  esac
done

shift $((OPTIND - 1))

gshuf -n $number $input_file > rand-${number}-${input_file}
