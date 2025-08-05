#!/bin/sh
USED=$(df / | awk 'NR==2 {print $5}' | tr -d '%')
if ["$USED" -gt 80]; then
	echo "ディスク容量が危険！使用率：$USED%"
else 
	echo "ディスク容量は正常です。使用率：$USED%"
fi
