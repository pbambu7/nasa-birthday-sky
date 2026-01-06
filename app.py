import streamlit as st
import requests
from datetime import date
from PIL import Image, ImageDraw
import io
import time
import random

from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, get_body
from astropy import units as u
from geopy.geocoders import Nominatim

# ---------------- PAGE CONFIG ----------------
st.set_page_config(
    page_title="What Was Above You When You Were Born?",
    layout="centered"
)

# ---------------- NASA API KEY ----------------
NASA_API_KEY = st.secrets.get("NASA_API_KEY")
if not NASA_API_KEY:
    st.error("NASA API key missing. Add it in Streamlit Secrets.")
    st.stop()

# ---------------- UI ----------------
st.title("🌌 What Was Above You When You Were Born?")
st.caption("Built with real NASA data")

st.markdown(
    "<small>🛰️ NASA API Status: <span style='color:orange'>Degraded</span></small>",
    unsafe_allow_html=True
)

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

    for _ in range(3):
        try:
            response = requests.get(url, params=params, timeout=8)
            response.raise_for_status()
            data = response.json()

            # Skip videos (YouTube/Vimeo)
            if data.get("media_type") != "image":
                return None

            return data

        except requests.exceptions.ReadTimeout:
            time.sleep(1)
        except requests.exceptions.RequestException:
            return None

    return None


@st.cache_data(show_spinner=False)
def get_visible_planets(birth_date, city):
    try:
        geolocator = Nominatim(
            user_agent="nasa-birthday-sky",
            timeout=5
        )
        location_data = geolocator.geocode(city)
    except Exception:
        return []

    if not location_data:
        return []

    location = EarthLocation(
        lat=location_data.latitude * u.deg,
        lon=location_data.longitude * u.deg
    )

    obs_time = Time(str(birth_date))
    frame = AltAz(obstime=obs_time, location=location)

    planets = ["mars", "venus", "jupiter"]
    visible = []

    for planet in planets:
        try:
            body = get_body(planet, obs_time).transform_to(frame)
            if body.alt > 0 * u.deg:
                visible.append(planet.capitalize())
        except Exception:
            continue

    return visible


def create_share_card(image_url, birth_date, city):
    try:
        img_bytes = requests.get(image_url, timeout=8).content
        img = Image.open(io.BytesIO(img_bytes)).convert("RGB")
    except Exception:
        return None

    img = img.resize((1080, 1080))
    draw = ImageDraw.Draw(img)

    draw.rectangle([(0, 820), (1080, 1080)], fill=(0, 0, 0))

    text = (
        "What was above me when I was born?\n"
        f"{birth_date}\n"
        f"{city}\n"
        "NASA Open Data"
    )

    draw.text((40, 860), text, fill="white")

    buf = io.BytesIO()
    img.save(buf, format="PNG")
    buf.seek(0)
    return buf


def generate_starfield():
    img = Image.new("RGB", (1080, 1080), "black")
    draw = ImageDraw.Draw(img)

    for _ in range(1200):
        x = random.randint(0, 1079)
        y = random.randint(0, 1079)
        brightness = random.randint(180, 255)
        draw.point((x, y), fill=(brightness, brightness, brightness))

    draw.text(
        (40, 900),
        "NASA APOD temporarily unavailable\nPlanet calculations still accurate",
        fill="white"
    )

    buf = io.BytesIO()
    img.save(buf, format="PNG")
    buf.seek(0)
    return buf

# ---------------- ACTION ----------------
if st.button("🚀 Run my birthday"):
    if not birth_city.strip():
        st.warning("Please enter a birth city.")
        st.stop()

    with st.spinner("Accessing NASA data..."):
        apod = fetch_apod(birth_date)

    if apod:
        st.subheader(apod.get("title", "NASA Astronomy Picture of the Day"))
        st.image(apod["url"], use_column_width=True)
        st.write(apod.get("explanation", ""))

        card = create_share_card(apod["url"], birth_date, birth_city)

    else:
        st.warning(
            "🚨 NASA APOD is temporarily unavailable.\n"
            "Planetary calculations are still accurate."
        )

        fallback_card = generate_starfield()
        st.image(fallback_card, caption="Generated starfield (fallback)")
        card = fallback_card

    planets = get_visible_planets(birth_date, birth_city)

    st.subheader("🪐 Planets visible above the horizon")
    if planets:
        st.write(", ".join(planets))
    else:
        st.write("No major planets were visible at that moment.")

    if card:
        st.download_button(
            "📸 Download IG Share Card",
            data=card,
            file_name="my_birth_universe.png",
            mime="image/png"
        )
