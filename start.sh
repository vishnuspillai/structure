#!/bin/bash

command -v python3 >/dev/null 2>&1 || { echo >&2 "Python3 required"; exit 1; }
command -v npm >/dev/null 2>&1 || { echo >&2 "Node.js required"; exit 1; }

echo "Setting up RAREMISS..."

# Create virtual environment
python3 -m venv venv
source venv/bin/activate

# Install backend dependencies
pip install --upgrade pip
pip install -r requirements.txt

# Start backend
echo "Starting backend..."
uvicorn src.api.main:app --host 0.0.0.0 --port 8000 &

# Move to frontend
cd ui

# Install frontend dependencies
npm install

# Start frontend
echo "Starting frontend..."
npm run dev &

echo "=========================================="
echo "          RAREMISS is running!            "
echo "=========================================="
echo "Frontend: http://localhost:5173"
echo "Backend:  http://localhost:8000"
echo "=========================================="
wait
