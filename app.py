# app.py — NASA Birthday Sky (Streamlit)
# ------------------------------------------------------------
# Personal astronomy using NASA Open APIs + Astropy
# Fixed for Astropy >= 5.x / Python 3.13
# ------------------------------------------------------------

import io
import os
import csv
import math
import time
import random
import textwrap
from pathlib import Path
from datetime import date, datetime, time as dtime

from reportlab.lib.pagesizes import LETTER
from reportlab.pdfgen import canvas
from reportlab.lib.units import inch

import streamlit as st
import requests
from PIL import Image, ImageDraw, ImageFont

import pytz
from timezonefinder import TimezoneFinder

from astropy.time import Time
from astropy.coordinates import (
    EarthLocation,
    AltAz,
    SkyCoord,
    get_body,
    get_sun,
    get_moon,
    solar_system_ephemeris,
)
from astropy.coordinates import moon_illumination
from astropy import units as u

from geopy.geocoders import Nominatim


# =========================
# CONFIG / PAGE
# =========================
st.set_page_config(page_title="NASA Birthday Sky", layout="wide")

NASA_API_KEY = st.secrets.get("NASA_API_KEY")
APP_BASE_URL = st.secrets.get("APP_BASE_URL", "").strip()

if not NASA_API_KEY:
    st.error(
        'NASA API key missing.\n\nAdd in Streamlit → App settings → Secrets:\n\nNASA_API_KEY = "YOUR_KEY"'
    )
    st.stop()

STAR_CATALOG_SOURCE = (
    "Gaia Collaboration et al., Gaia Data Release 3 (DR3), "
    "European Space Agency (ESA), 2022. "
    "https://gea.esac.esa.int/archive/"
)

DATA_DIR = Path("data")
DATA_DIR.mkdir(exist_ok=True)
DEFAULT_CSV_PATH = DATA_DIR / "bright_stars.csv"

st.title("🌌 NASA Birthday Sky")
st.caption("Personal astronomy using NASA open data + real sky calculations (timezone-correct).")


# =========================
# EMBEDDED BRIGHT STARS
# =========================
BRIGHT_STARS = {
    "Sirius": (101.2875, -16.7161, -1.46),
    "Canopus": (95.9879, -52.6957, -0.74),
    "Arcturus": (213.9153, 19.1824, -0.05),
    "Vega": (279.2347, 38.7837, 0.03),
    "Capella": (79.1723, 45.9979, 0.08),
    "Rigel": (78.6345, -8.2016, 0.12),
    "Procyon": (114.8255, 5.2250, 0.38),
    "Betelgeuse": (88.7929, 7.4071, 0.42),
    "Altair": (297.6958, 8.8683, 0.77),
}


# =========================
# STAR CSV UTILITIES
# =========================
def load_star_catalog(path: Path, max_rows=5000):
    stars = []
    try:
        with open(path, newline="") as f:
            reader = csv.DictReader(f)
            for i, row in enumerate(reader):
                if i >= max_rows:
                    break
                stars.append(
                    {
                        "name": row.get("name", f"Star {i+1}"),
                        "ra": float(row["ra"]),
                        "dec": float(row["dec"]),
                        "mag": float(row.get("mag", 0.0)),
                    }
                )
    except Exception:
        return []
    return stars


# =========================
# INPUTS
# =========================
col1, col2 = st.columns(2)
with col1:
    birth_date = st.date_input("Your birth date", value=date(2000, 1, 1))
with col2:
    birth_city = st.text_input("Birth city", value="Salt Lake City")

moment = st.selectbox(
    "Moment to calculate (local time)",
    ["Midnight (00:00)", "Noon (12:00)", "Exact time"],
)

birth_time = dtime(0, 0)
if moment == "Noon (12:00)":
    birth_time = dtime(12, 0)
elif moment == "Exact time":
    birth_time = st.time_input("Exact birth time (HH:MM)", value=dtime(0, 0))

educator_mode = st.toggle("🧪 Educator Mode (guided lesson + explain the science)")

run_clicked = st.button("🚀 Run my birthday")
if not run_clicked:
    st.stop()


# =========================
# GEO + TIMEZONE
# =========================
@st.cache_data(show_spinner=False)
def geocode_city(city):
    geolocator = Nominatim(user_agent="nasa-birthday-sky")
    loc = geolocator.geocode(city)
    if not loc:
        return None
    return float(loc.latitude), float(loc.longitude), loc.address


geo = geocode_city(birth_city)
if not geo:
    st.error("Could not geocode that city.")
    st.stop()

lat, lon, address = geo
location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg)

tf = TimezoneFinder()
tz_name = tf.timezone_at(lat=lat, lng=lon)
tz = pytz.timezone(tz_name)

local_dt = tz.localize(datetime.combine(birth_date, birth_time))
utc_dt = local_dt.astimezone(pytz.utc)
obs_time = Time(utc_dt)


# =========================
# PLANETS
# =========================
def compute_visible_planets(location, obs_time):
    frame = AltAz(obstime=obs_time, location=location)
    planets = []
    for p in ["venus", "mars", "jupiter"]:
        body = get_body(p, obs_time).transform_to(frame)
        if body.alt.deg > 0:
            planets.append(
                {"name": p.capitalize(), "alt": body.alt.deg, "az": body.az.deg}
            )
    return planets


# =========================
# 🌕 FIXED MOON PHASE (ASTROPY-SAFE)
# =========================
def compute_moon_phase(obs_time):
    illum = moon_illumination(obs_time)

    with solar_system_ephemeris.set("builtin"):
        moon = get_moon(obs_time)
        sun = get_sun(obs_time)

    elongation = moon.separation(sun).deg
    waxing = elongation < 180
    pct = illum * 100

    if pct < 1:
        name = "New Moon"
    elif pct < 49:
        name = "Waxing Crescent" if waxing else "Waning Crescent"
    elif 49 <= pct <= 51:
        name = "First Quarter" if waxing else "Last Quarter"
    elif pct < 99:
        name = "Waxing Gibbous" if waxing else "Waning Gibbous"
    else:
        name = "Full Moon"

    return float(illum), name, waxing


# =========================
# STARS
# =========================
def compute_visible_stars(obs_time, location):
    visible = []
    for name, (ra, dec, mag) in BRIGHT_STARS.items():
        sc = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)
        altaz = sc.transform_to(AltAz(obstime=obs_time, location=location))
        if altaz.alt.deg > 0:
            visible.append(
                {
                    "name": name,
                    "alt": altaz.alt.deg,
                    "az": altaz.az.deg,
                    "mag": mag,
                }
            )
    return visible


# =========================
# RUN COMPUTATION
# =========================
planets = compute_visible_planets(location, obs_time)
moon_illum, moon_phase_name, waxing = compute_moon_phase(obs_time)
stars = compute_visible_stars(obs_time, location)

moon_pct = int(round(moon_illum * 100))

st.subheader("🪐 Planets Visible")
if planets:
    for p in planets:
        st.write(f"**{p['name']}** — Alt {p['alt']:.1f}°, Az {p['az']:.1f}°")
else:
    st.write("No major planets visible.")

st.subheader("🌕 Moon")
st.write(f"**{moon_phase_name}** — {moon_pct}% illuminated")

st.subheader("⭐ Bright Stars")
st.write(f"{len(stars)} stars above the horizon")

if educator_mode:
    st.info(
        "Educator Mode: Altitude > 0° means visible above the horizon. "
        "Moon illumination is computed using Sun–Moon phase angle (Astropy)."
    )

st.caption(
    f"📍 {address} • Timezone: {tz_name} • UTC used: {utc_dt.isoformat()}"
)
