#!/bin/bash
echo "Stopping RAREMISS..."

pkill -f uvicorn
pkill -f node

echo "Stopped."
