🌌 NASA Birthday Sky
What did the sky look like on the day you were born?
NASA Birthday Sky is a simple interactive app that lets anyone enter their birthdate and location to see which planets (Mars, Venus, Jupiter) were visible in the sky at that moment — using real NASA data.
Built to make space science personal, visual, and shareable.
🚀 Features
🌍 Location-based sky accuracy
📅 Birthdate input (run your birthday)
🪐 Planet visibility (Mars, Venus, Jupiter)
🖼️ Shareable, IG-ready sky image
🔬 Powered by official NASA APIs
🧠 Why this project
I was inspired by seeing how AI and open science can unlock discovery — not just for researchers, but for everyday people.
This project explores how public NASA data can be turned into an experience that’s:
Educational
Personal
Visually engaging
Built with a focus on clarity, accuracy, and execution.
🛠️ Tech Stack
Python
Streamlit
NASA APIs
Astropy
Geopy
Pillow
🔐 API Key Setup
This app uses a NASA API key.
If running locally:
Get a free key from: https://api.nasa.gov
Set it as an environment variable:
export NASA_API_KEY=your_api_key_here
Streamlit Cloud users should add the key under Secrets:
NASA_API_KEY = your_api_key_here
▶️ Run Locally
pip install -r requirements.txt
streamlit run app.py
🌍 Live Demo
👉 https://nasa-birthday-sky.streamlit.app (once deployed)
📌 Notes
This project focuses on planets, not asteroids (yet)
Designed for public sharing and accessibility
Educational use — not intended for navigation or observation planning
🙌 Credits
Data provided by NASA Open APIs
Built by Patience Bambu
