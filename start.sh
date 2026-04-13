#!/bin/bash
set -e

command -v python3 >/dev/null 2>&1 || { echo >&2 "Python3 required"; exit 1; }
command -v npm >/dev/null 2>&1 || { echo >&2 "Node.js required"; exit 1; }

if lsof -i:8000 >/dev/null; then
  echo "Port 8000 already in use. Please free it."
  exit 1
fi

if lsof -i:5173 >/dev/null; then
  echo "Port 5173 already in use. Please free it."
  exit 1
fi

echo "Setting up RAREMISS..."

if [ -d "venv" ]; then
  echo "Removing old virtual environment..."
  rm -rf venv
fi

# Create virtual environment
python3 -m venv venv
source venv/bin/activate

# Install backend dependencies
pip install --upgrade pip
pip install -r requirements.txt

# Start backend
echo "Starting backend..."
uvicorn src.api.main:app --host 0.0.0.0 --port 8000 > backend.log 2>&1 &

# Move to frontend
cd ui

# Install frontend dependencies
npm install

# Start frontend
echo "Starting frontend..."
npm run dev > ../frontend.log 2>&1 &
cd ..

echo "===================================="
echo "RAREMISS is running"
echo "Frontend: http://localhost:5173"
echo "Backend:  http://localhost:8000"
echo "Logs:"
echo "  backend.log"
echo "  frontend.log"
echo "===================================="
wait
