#!/bin/bash

# This script generates git commit logs between two tags
# Usage: ./tag-diff.sh <tag1> <tag2>
# 
# To create a changelog section from the output:
# 1. Run this script with two version tags (e.g., ./tag-diff.sh 0.1.8 0.1.9)
# 2. Copy the output and ask Claude to "write it as a keepachangelog section in md"
# 3. Claude will convert the git commits into proper Keep a Changelog format
#    with Added/Changed/Fixed sections following keepachangelog.com standards

if [ $# -ne 2 ]; then
    echo "Usage: $0 <tag1> <tag2>"
    echo "Example: $0 0.1.8 0.1.9"
    exit 1
fi

TAG1=$1
TAG2=$2

echo "Fetching all tags..."
git fetch --tags

echo "Collecting commits between $TAG1 and $TAG2..."
git log $TAG1..$TAG2

