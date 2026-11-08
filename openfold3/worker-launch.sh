#!/bin/zsh
set -e
cd /Users/yili/Documents/ChatGPT/autodock
exec /Users/yili/Documents/ChatGPT/autodock/.env/bin/python3.11 /Users/yili/Documents/ChatGPT/autodock/openfold3/worker.py --host 100.64.0.22 --port 8099 --token-file /Users/yili/.config/marine-agent/openfold-worker-token
