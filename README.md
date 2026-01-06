🌌 NASA Birthday Sky
What did the sky look like on the day you were born?
NASA Birthday Sky is an interactive web app that lets anyone enter their birthdate and location to discover which major planets — Mars, Venus, and Jupiter — were visible in the sky at that moment, using real NASA open data and astronomy calculations.
The goal is simple: make space science personal, visual, and shareable.
🚀 Features
🌍 Location-based sky accuracy using real geographic coordinates
📅 Birthdate input — run your birthday in seconds
🪐 Planet visibility detection (Mars, Venus, Jupiter)
🖼️ Shareable, IG-ready sky image
🔬 Powered by official NASA Open APIs
🧯 Resilient design — app remains usable even during NASA API outages
🧠 Why This Project
NASA’s open data is incredibly powerful, but often inaccessible to people without a technical background.
This project explores how open science, combined with thoughtful engineering and design, can turn raw scientific data into an experience that is:
Educational
Personal
Visually engaging
Built with a focus on clarity, accuracy, resilience, and execution, NASA Birthday Sky is designed for public engagement — not just research labs.
🛠️ Tech Stack
Python
Streamlit
NASA Open APIs (APOD)
Astropy — astronomical calculations
Geopy — location → coordinates
Pillow — image generation & share cards
🔐 API Key Setup
This app requires a NASA API key.
Get a free key:
👉 https://api.nasa.gov
Run locally:
export NASA_API_KEY=your_api_key_here
Streamlit Cloud:
Add the key under App Settings → Secrets:
NASA_API_KEY = "your_api_key_here"

▶️ Run Locally
pip install -r requirements.txt
streamlit run app.py
