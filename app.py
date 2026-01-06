import streamlit as st
import requests
from datetime import date
from PIL import Image, ImageDraw
import io
import os

from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, get_body
from astropy import units as u
from geopy.geocoders import Nominatim

# ---------------- CONFIG ----------------
st.set_page_config(
    page_title="What Was Above You When You Were Born?",
    layout="centered"
)

# 🔐 Secure API key handling
NASA_API_KEY = st.secrets.get("NASA_API_KEY") or os.getenv("NASA_API_KEY")

if not NASA_API_KEY:
    st.error("NASA API key not found. Please set NASA_API_KEY in Streamlit Secrets.")
    st.stop()

# ---------------- UI ----------------
st.title("🌌 What Was Above You When You Were Born?")
st.caption("Built with real NASA data")

birth_date = st.date_input("Your birth date", value=date(2000, 1, 1))
birth_city = st.text_input("Birth city", placeholder="Salt Lake City")

# ---------------- FUNCTIONS ----------------
@st.cache_data(show_spinner=False)
def fetch_apod(birth_date):
    url = "https://api.nasa.gov/planetary/apod"
    params = {
        "date": birth_date.isoformat(),
        "api_key": NASA_API_KEY
    }
    response = requests.get(url, params=params, timeout=20)
    response.raise_for_status()
    return response.json()

@st.cache_data(show_spinner=False)
def get_visible_planets(birth_date, city):
    geolocator = Nominatim(user_agent="nasa_birthday_sky_app")
    location_data = geolocator.geocode(city)

    if not location_data:
        return []

    location = EarthLocation(
        lat=location_data.latitude * u.deg,
        lon=location_data.longitude * u.deg
    )

    time = Time(str(birth_date))
    frame = AltAz(obstime=time, location=location)

    planets = ["mars", "venus", "jupiter"]
    visible = []

    for planet in planets:
        body = get_body(planet, time).transform_to(frame)
        if body.alt > 0 * u.deg:
            visible.append(planet.capitalize())

    return visible

def create_share_card(image_url, birth_date, city):
    img_bytes = requests.get(image_url, timeout=20).content
    img = Image.open(io.BytesIO(img_bytes)).convert("RGB")
    img = img.resize((1080, 1080))

    draw = ImageDraw.Draw(img)

    overlay_height = 260
    draw.rectangle(
        [(0, 1080 - overlay_height), (1080, 1080)],
        fill=(0, 0, 0)
    )

    text = (
        "What was above me when I was born?\n"
        f"{birth_date}\n"
        f"{city}\n"
        "NASA Open Data"
    )

    draw.text((50, 820), text, fill="white")

    buffer = io.BytesIO()
    img.save(buffer, format="PNG")
    buffer.seek(0)
    return buffer

# ---------------- ACTION ----------------
if st.button("🚀 Run my birthday"):
    if not birth_city.strip():
        st.warning("Please enter a birth city.")
    else:
        try:
            with st.spinner("Accessing NASA data..."):
                apod = fetch_apod(birth_date)

            st.subheader(apod.get("title", "NASA Astronomy Picture of the Day"))

            if apod.get("media_type") == "image":
                st.image(apod["url"], use_column_width=True)
            else:
                st.write("Media type is not an image for this date.")

            st.write(apod.get("explanation", ""))

            planets = get_visible_planets(birth_date, birth_city)

            st.subheader("🪐 Planets visible above the horizon")
            if planets:
                st.write(", ".join(planets))
            else:
                st.write("No major planets were visible at that moment.")

            card = create_share_card(apod["url"], birth_date, birth_city)

            st.download_button(
                "📸 Download IG Share Card",
                data=card,
                file_name="my_birth_universe.png",
                mime="image/png"
            )

        except Exception as e:
            st.error("Something went wrong while fetching NASA data.")
            st.exception(e)

# ---------------- FOOTER ----------------
st.caption(
    "Educational use only · Data from NASA Open APIs · "
    "Not intended for navigation or astronomical planning\n\n"
    "Built by Patience Bambu 🚀"
)
