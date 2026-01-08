# app.py — NASA Birthday Sky (Streamlit)
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
    get_body,
    SkyCoord,
    BarycentricTrueEcliptic,
)
from astropy import units as u

from geopy.geocoders import Nominatim


# =========================
# CONFIG / PAGE
# =========================
st.set_page_config(page_title="NASA Birthday Sky", layout="wide")

st.markdown(
    """
<style>
.block-container { padding-top: 2rem; padding-bottom: 3rem; max-width: 980px; }
.stButton>button { width:100%; padding:14px; border-radius:14px; font-size:16px; }
small.muted { color: rgba(255,255,255,0.65); }
hr { margin: 1.25rem 0; }

.card {
  border: 1px solid rgba(255,255,255,0.10);
  background: rgba(255,255,255,0.04);
  border-radius: 16px;
  padding: 16px;
}
.card h4 { margin: 0 0 6px 0; }
.card p { margin: 0; color: rgba(255,255,255,0.82); line-height: 1.35; }
.badge {
  display:inline-block;
  padding: 4px 10px;
  border-radius:999px;
  background: rgba(255,255,255,0.08);
  border: 1px solid rgba(255,255,255,0.14);
  font-size: 12px;
  margin-right: 8px;
}

.kpi { display:flex; gap:12px; flex-wrap:wrap; }
.kpi .pill {
  border: 1px solid rgba(255,255,255,0.12);
  background: rgba(255,255,255,0.05);
  border-radius: 999px;
  padding: 8px 12px;
  font-size: 13px;
  color: rgba(255,255,255,0.86);
}
</style>
""",
    unsafe_allow_html=True,
)

NASA_API_KEY = st.secrets.get("NASA_API_KEY")
APP_BASE_URL = st.secrets.get("APP_BASE_URL", "").strip()

STAR_CATALOG_SOURCE = st.secrets.get(
    "STAR_CATALOG_SOURCE",
    "Embedded bright-star subset (RA/Dec). For production-grade catalogs, use Hipparcos / Tycho-2 or Gaia DR3.",
)

DATA_DIR = Path("data")
DATA_DIR.mkdir(exist_ok=True)
DEFAULT_CSV_PATH = DATA_DIR / "bright_stars.csv"

if not NASA_API_KEY:
    st.error('NASA API key missing. Add it in Streamlit → App settings → Secrets:\n\nNASA_API_KEY = "YOUR_KEY"')
    st.stop()

st.title("🌌 NASA Birthday Sky")
st.caption("Personal astronomy using NASA open data + real sky calculations (timezone-correct).")


# =========================
# EMBEDDED BRIGHT STARS (compact subset)
# =========================
# RA/Dec in degrees (ICRS).
# NOTE: mag=0.0 here; CSV mode supports mag if provided.
BRIGHT_STARS = {
    # Orion
    "Betelgeuse": (88.7929, 7.4071, 0.42),
    "Rigel": (78.6345, -8.2016, 0.13),
    "Bellatrix": (81.2828, 6.3497, 1.64),
    "Saiph": (86.9391, -9.6696, 2.06),
    "Alnitak": (85.1897, -1.9426, 1.74),
    "Alnilam": (84.0534, -1.2019, 1.69),
    "Mintaka": (83.0017, -0.2991, 2.25),

    # Ursa Major
    "Dubhe": (165.9322, 61.7510, 1.79),
    "Merak": (165.4603, 56.3824, 2.34),
    "Phecda": (178.4577, 53.6948, 2.41),
    "Megrez": (183.8565, 57.0326, 3.32),
    "Alioth": (193.5073, 55.9598, 1.76),
    "Mizar": (200.9814, 54.9254, 2.23),
    "Alkaid": (206.8856, 49.3133, 1.85),

    # Leo
    "Regulus": (152.0929, 11.9672, 1.35),
    "Denebola": (177.2649, 14.5721, 2.14),
    "Algieba": (154.9931, 19.8415, 2.01),

    # Summer triangle
    "Vega": (279.2347, 38.7837, 0.03),
    "Deneb": (310.3579, 45.2803, 1.25),
    "Altair": (297.6958, 8.8683, 0.77),

    # Others
    "Sirius": (101.2875, -16.7161, -1.46),
    "Canopus": (95.9879, -52.6957, -0.74),
    "Arcturus": (213.9153, 19.1824, -0.05),
    "Spica": (201.2983, -11.1613, 0.98),
    "Antares": (247.3519, -26.4320, 1.06),
    "Capella": (79.1723, 45.9979, 0.08),
}

# Constellation line segments defined by star names in BRIGHT_STARS
CONSTELLATION_LINES = {
    "Orion": [
        ("Betelgeuse", "Bellatrix"),
        ("Bellatrix", "Mintaka"),
        ("Mintaka", "Alnilam"),
        ("Alnilam", "Alnitak"),
        ("Alnitak", "Saiph"),
        ("Saiph", "Rigel"),
        ("Rigel", "Betelgeuse"),
    ],
    "Ursa Major": [
        ("Dubhe", "Merak"),
        ("Merak", "Phecda"),
        ("Phecda", "Megrez"),
        ("Megrez", "Alioth"),
        ("Alioth", "Mizar"),
        ("Mizar", "Alkaid"),
    ],
    "Leo": [
        ("Regulus", "Algieba"),
        ("Algieba", "Denebola"),
    ],
    "Summer Triangle": [
        ("Vega", "Deneb"),
        ("Deneb", "Altair"),
        ("Altair", "Vega"),
    ],
}

# =========================
# STAR CSV UTILITIES
# =========================
def write_default_bright_stars_csv(path: Path) -> None:
    """Create/overwrite data/bright_stars.csv from embedded BRIGHT_STARS."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["name", "ra", "dec", "mag"])
        for name, (ra, dec, mag) in sorted(BRIGHT_STARS.items(), key=lambda x: x[0].lower()):
            w.writerow([name, f"{ra:.6f}", f"{dec:.6f}", f"{mag:.2f}"])


def load_star_catalog(path: Path, max_rows=5000):
    """
    CSV catalog loader.
    Expected columns: name, ra, dec, mag (mag optional)
    RA/DEC in degrees.
    """
    stars = []
    try:
        with open(path, newline="") as f:
            reader = csv.DictReader(f)
            for i, row in enumerate(reader):
                if i >= max_rows:
                    break
                ra = float(row["ra"])
                dec = float(row["dec"])
                name = (row.get("name") or "").strip() or f"Star {i+1}"
                mag_raw = row.get("mag", "")
                mag = float(mag_raw) if str(mag_raw).strip() not in ("", "None", "nan") else 0.0
                stars.append({"name": name, "ra": ra, "dec": dec, "mag": mag})
    except Exception:
        return []
    return stars


# =========================
# QUERY PARAMS (SHARE LINKS)
# =========================
def get_query_params():
    try:
        return dict(st.query_params)
    except Exception:
        try:
            qp = st.experimental_get_query_params()
            return {k: v[0] if isinstance(v, list) and v else v for k, v in qp.items()}
        except Exception:
            return {}

qp = get_query_params()

def qp_get(key, default=""):
    v = qp.get(key, default)
    if isinstance(v, list):
        return v[0] if v else default
    return v if v is not None else default

def parse_iso_date(s, fallback):
    try:
        return date.fromisoformat(str(s))
    except Exception:
        return fallback

def parse_hhmm(s, fallback_time=dtime(0, 0)):
    try:
        parts = str(s).split(":")
        hh = max(0, min(23, int(parts[0])))
        mm = max(0, min(59, int(parts[1]) if len(parts) > 1 else 0))
        return dtime(hh, mm)
    except Exception:
        return fallback_time

default_birth_date = parse_iso_date(qp_get("date", "2000-01-01"), date(2000, 1, 1))
default_city = qp_get("city", "")
default_mode = qp_get("mode", "midnight").lower()
default_time = parse_hhmm(qp_get("time", "00:00"), dtime(0, 0))
default_autorun = str(qp_get("autorun", "0")) == "1"


# =========================
# SIDEBAR: CATALOG + RENDER SETTINGS
# =========================
with st.sidebar:
    st.header("⚙️ Settings")

    st.subheader("Star Catalog")
    catalog_mode = st.radio(
        "Catalog source",
        options=["Embedded (fast)", "CSV (data/bright_stars.csv)"],
        index=0,
    )

    st.caption("If you choose CSV mode, the file must exist with columns: name, ra, dec, mag (optional).")

    col_a, col_b = st.columns(2)
    with col_a:
        if st.button("📝 Create data/bright_stars.csv", use_container_width=True):
            write_default_bright_stars_csv(DEFAULT_CSV_PATH)
            st.success(f"Created: {DEFAULT_CSV_PATH.as_posix()}")
    with col_b:
        csv_exists = DEFAULT_CSV_PATH.exists()
        st.metric("CSV status", "Found ✅" if csv_exists else "Missing ❌")

    st.divider()

    st.subheader("Rendering")
    show_constellations = st.toggle("Draw constellation lines", value=True)
    label_stars = st.toggle("Label top stars", value=True)
    max_render_stars = st.slider("Max stars to render", 50, 1200, 300, step=50)
    max_compute_stars = st.slider("Max stars to compute (speed)", 100, 5000, 1500, step=100)

    st.caption("LOD = ‘Level of Detail’: render fewer points for speed + clarity.")

    st.divider()
    st.subheader("PDF")
    include_apod_excerpt = st.toggle("Include APOD explanation excerpt in PDF", value=False)


# =========================
# UI INPUTS
# =========================
top1, top2 = st.columns(2)
with top1:
    birth_date = st.date_input("Your birth date", value=default_birth_date)
with top2:
    birth_city = st.text_input("Birth city", value=default_city, placeholder="Salt Lake City, UT")

mode_labels = {"midnight": "Midnight (00:00)", "noon": "Noon (12:00)", "exact": "Exact time"}
mode_options = ["Midnight (00:00)", "Noon (12:00)", "Exact time"]
default_mode_label = mode_labels.get(default_mode, "Midnight (00:00)")

moment = st.selectbox("Moment to calculate (local time)", mode_options, index=mode_options.index(default_mode_label))

birth_time = dtime(0, 0)
if moment.startswith("Midnight"):
    birth_time = dtime(0, 0)
elif moment.startswith("Noon"):
    birth_time = dtime(12, 0)
else:
    birth_time = st.time_input("Exact birth time (HH:MM) — optional", value=default_time)

educator_mode = st.toggle("🧪 Educator Mode (guided lesson + explain the science)", value=False)

st.divider()


# =========================
# CACHE HELPERS
# =========================
@st.cache_data(show_spinner=False, ttl=60 * 60 * 24)
def geocode_city(city: str):
    city = city.strip()
    if not city:
        return None
    try:
        geolocator = Nominatim(user_agent="nasa-birthday-sky", timeout=6)
        loc = geolocator.geocode(city)
        if not loc:
            return None
        return float(loc.latitude), float(loc.longitude), str(loc.address)
    except Exception:
        return None


def local_dt_to_utc(birth_date: date, birth_time: dtime, lat: float, lon: float):
    tf = TimezoneFinder()
    tz_name = tf.timezone_at(lat=lat, lng=lon)
    if not tz_name:
        return None, None

    tz = pytz.timezone(tz_name)
    local_dt = tz.localize(datetime.combine(birth_date, birth_time))
    utc_dt = local_dt.astimezone(pytz.utc)
    return tz_name, utc_dt


@st.cache_data(show_spinner=False, ttl=60 * 60)
def fetch_apod(date_iso: str, api_key: str):
    url = "https://api.nasa.gov/planetary/apod"
    params = {"date": date_iso, "api_key": api_key}
    for attempt in range(4):
        try:
            r = requests.get(url, params=params, timeout=(4, 10))
            r.raise_for_status()
            return r.json()
        except requests.exceptions.Timeout:
            time.sleep(0.6 * (2 ** attempt))
        except Exception:
            return None
    return None


# =========================
# ASTRONOMY: PLANETS + MOON + STARS
# =========================
def compute_visible_planets(location: EarthLocation, obs_time: Time):
    frame = AltAz(obstime=obs_time, location=location)
    planet_names = ["venus", "mars", "jupiter"]
    results = []
    for p in planet_names:
        try:
            body = get_body(p, obs_time).transform_to(frame)
            alt = body.alt.to(u.deg).value
            az = body.az.to(u.deg).value
            if alt > 0:
                results.append({"name": p.capitalize(), "alt": float(alt), "az": float(az)})
        except Exception:
            continue
    results.sort(key=lambda x: x["alt"], reverse=True)
    return results


def compute_moon_phase(obs_time: Time):
    moon = get_body("moon", obs_time)
    sun = get_body("sun", obs_time)

    D = moon.separation(sun).to(u.rad).value
    illum = (1 - math.cos(D)) / 2

    moon_ecl = moon.transform_to(BarycentricTrueEcliptic(obstime=obs_time))
    sun_ecl = sun.transform_to(BarycentricTrueEcliptic(obstime=obs_time))
    dlon = (moon_ecl.lon - sun_ecl.lon).wrap_at(360 * u.deg).to(u.deg).value
    waxing = dlon > 0

    pct = illum * 100
    if pct < 2:
        name = "New Moon"
    elif pct < 48:
        name = "Waxing Crescent" if waxing else "Waning Crescent"
    elif pct < 52:
        name = "First Quarter" if waxing else "Last Quarter"
    elif pct < 98:
        name = "Waxing Gibbous" if waxing else "Waning Gibbous"
    else:
        name = "Full Moon"
    return float(illum), name, waxing


def project_star_to_altaz(star_ra_deg, star_dec_deg, obs_time: Time, location: EarthLocation):
    sc = SkyCoord(ra=star_ra_deg * u.deg, dec=star_dec_deg * u.deg, frame="icrs")
    altaz = sc.transform_to(AltAz(obstime=obs_time, location=location))
    return altaz.alt.to(u.deg).value, altaz.az.to(u.deg).value


@st.cache_data(show_spinner=False, ttl=60 * 60)
def compute_visible_stars_embedded(obs_time_iso: str, lat: float, lon: float):
    obs_time = Time(obs_time_iso)
    location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg)
    visible = []
    for name, (ra, dec, mag) in BRIGHT_STARS.items():
        alt, az = project_star_to_altaz(ra, dec, obs_time, location)
        if alt > 0:
            visible.append({"name": name, "ra": ra, "dec": dec, "alt": float(alt), "az": float(az), "mag": float(mag)})
    visible.sort(key=lambda s: s["alt"], reverse=True)
    return visible


@st.cache_data(show_spinner=False, ttl=60 * 60)
def compute_visible_stars_csv(obs_time_iso: str, lat: float, lon: float, csv_path_str: str, max_rows=1500):
    csv_path = Path(csv_path_str)
    catalog = load_star_catalog(csv_path, max_rows=max_rows)
    if not catalog:
        return []

    obs_time = Time(obs_time_iso)
    location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg)

    visible = []
    for star in catalog:
        try:
            alt, az = project_star_to_altaz(star["ra"], star["dec"], obs_time, location)
            if alt > 0:
                visible.append({
                    "name": star["name"],
                    "ra": float(star["ra"]),
                    "dec": float(star["dec"]),
                    "mag": float(star.get("mag", 0.0)),
                    "alt": float(alt),
                    "az": float(az),
                })
        except Exception:
            continue
    visible.sort(key=lambda s: s["alt"], reverse=True)
    return visible


# =========================
# PDF: auto citations + pagination
# =========================
def _wrap_lines(text: str, width=95):
    return textwrap.wrap(text, width=width)

def _pdf_new_page(c):
    c.showPage()
    c.setFont("Helvetica", 11)

def _pdf_ensure_space(c, y, need=24):
    # if not enough space, start new page
    if y < 0.9 * inch + need:
        _pdf_new_page(c)
        return LETTER[1] - 1 * inch
    return y

def _pdf_draw_wrapped(c, x, y, text, width=95, line_height=14):
    for line in _wrap_lines(text, width=width):
        y = _pdf_ensure_space(c, y, need=line_height + 8)
        c.drawString(x, y, line)
        y -= line_height
    return y

class CitationManager:
    def __init__(self):
        self._order = []
        self._map = {}

    def cite(self, key: str, entry: str) -> str:
        """Return citation marker like [1], registering entry if new."""
        if key not in self._map:
            self._order.append((key, entry))
            self._map[key] = len(self._order)
        return f"[{self._map[key]}]"

    def entries(self):
        # returns list of (n, entry)
        out = []
        for key, entry in self._order:
            out.append((self._map[key], entry))
        out.sort(key=lambda x: x[0])
        return out


def generate_science_report_buffer(
    birth_date_str,
    city,
    lat,
    lon,
    timezone,
    utc_dt_str,
    moment_label,
    planets,
    moon_text,
    repro_link,
    apod_status_text,
    apod_title,
    apod_expl,
    stars_used,
    star_source_text,
    catalog_mode_used,
    include_apod_excerpt=False,
):
    buf = io.BytesIO()
    c = canvas.Canvas(buf, pagesize=LETTER)
    width, height = LETTER

    cm = CitationManager()

    # Define citation entries (stable, not dependent on live web)
    nasa_api = cm.cite("nasa_api", "NASA Open APIs (APOD): https://api.nasa.gov")
    astropy = cm.cite("astropy", "Astropy project (coordinates, transforms, ephemerides): https://www.astropy.org")
    nominatim = cm.cite("nominatim", "OpenStreetMap Nominatim geocoding (via Geopy).")
    tzfinder = cm.cite("timezonefinder", "Timezone lookup: timezonefinder + pytz.")
    # Catalog references (generic, safe)
    gaia = cm.cite("gaia_dr3", "Gaia DR3 (ESA): for production-grade star catalogs and magnitudes.")
    hip = cm.cite("hipparcos", "Hipparcos / Tycho catalogs (ESA): classic bright-star references for outreach/education.")

    y = height - 1 * inch
    c.setFont("Helvetica-Bold", 16)
    c.drawString(1 * inch, y, "NASA Birthday Sky — Science Report")

    y -= 0.45 * inch
    c.setFont("Helvetica", 11)
    c.drawString(1 * inch, y, f"Generated: {datetime.utcnow().isoformat()} UTC")

    # Inputs
    y -= 0.5 * inch
    c.setFont("Helvetica-Bold", 12)
    c.drawString(1 * inch, y, "Inputs")

    y -= 0.3 * inch
    c.setFont("Helvetica", 11)
    y = _pdf_draw_wrapped(c, 1 * inch, y, f"Birth date: {birth_date_str}")
    y -= 4
    y = _pdf_draw_wrapped(c, 1 * inch, y, f"City: {city} {nominatim}")
    y -= 4
    y = _pdf_draw_wrapped(c, 1 * inch, y, f"Lat/Lon: {lat:.4f}, {lon:.4f}")
    y -= 4
    y = _pdf_draw_wrapped(c, 1 * inch, y, f"Moment: {moment_label} (local time) → Timezone: {timezone} {tzfinder}")
    y -= 4
    y = _pdf_draw_wrapped(c, 1 * inch, y, f"UTC time used: {utc_dt_str}")

    # Results
    y -= 0.35 * inch
    c.setFont("Helvetica-Bold", 12)
    c.drawString(1 * inch, y, "Results")

    y -= 0.3 * inch
    c.setFont("Helvetica", 11)
    y = _pdf_draw_wrapped(c, 1 * inch, y, f"Moon: {moon_text} {astropy}")

    if planets:
        for p in planets:
            y -= 6
            y = _pdf_draw_wrapped(c, 1 * inch, y, f"{p['name']} — altitude {p['alt']:.2f}°, azimuth {p['az']:.2f}° {astropy}")
    else:
        y -= 6
        y = _pdf_draw_wrapped(c, 1 * inch, y, f"No major planets (Mars/Venus/Jupiter list) above the horizon at this moment. {astropy}")

    # Stars summary
    y -= 0.35 * inch
    c.setFont("Helvetica-Bold", 12)
    c.drawString(1 * inch, y, "Stars Used (Sky Overlay)")

    y -= 0.3 * inch
    c.setFont("Helvetica", 11)

    # cite catalog type
    if catalog_mode_used.startswith("CSV"):
        cat_cite = gaia  # “production-grade” pointer
    else:
        cat_cite = hip

    y = _pdf_draw_wrapped(
        c, 1 * inch, y,
        f"Catalog mode: {catalog_mode_used}. Source note: {star_source_text} {cat_cite}"
    )
    y -= 8

    if stars_used:
        y = _pdf_draw_wrapped(
            c, 1 * inch, y,
            f"Visible stars computed: {len(stars_used)} (sample below; alt>0°). {astropy}"
        )
        sample = stars_used[:14]
        for s in sample:
            y -= 4
            y = _pdf_draw_wrapped(
                c, 1 * inch, y,
                f"- {s['name']} — alt {s['alt']:.1f}°, az {s['az']:.1f}° (RA {s['ra']:.2f}°, Dec {s['dec']:.2f}°, mag {s.get('mag', 0.0):.2f})"
            )
    else:
        y = _pdf_draw_wrapped(c, 1 * inch, y, "No stars were computed/visible in the selected mode.")

    # APOD status
    y -= 0.35 * inch
    c.setFont("Helvetica-Bold", 12)
    c.drawString(1 * inch, y, "NASA APOD Status")

    y -= 0.3 * inch
    c.setFont("Helvetica", 11)
    y = _pdf_draw_wrapped(c, 1 * inch, y, f"{apod_status_text} {nasa_api}")

    if include_apod_excerpt and apod_expl:
        y -= 10
        c.setFont("Helvetica-Bold", 11)
        y = _pdf_draw_wrapped(c, 1 * inch, y, "APOD explanation excerpt:")
        c.setFont("Helvetica", 10)
        excerpt = apod_expl[:700].strip()
        y = _pdf_draw_wrapped(c, 1 * inch, y, excerpt, width=105, line_height=12)

    # Repro link
    y -= 0.35 * inch
    c.setFont("Helvetica-Bold", 12)
    c.drawString(1 * inch, y, "Reproducibility Link")

    y -= 0.3 * inch
    c.setFont("Helvetica", 11)
    y = _pdf_draw_wrapped(c, 1 * inch, y, repro_link)

    # Citations
    y -= 0.35 * inch
    c.setFont("Helvetica-Bold", 12)
    c.drawString(1 * inch, y, "References (Auto-numbered)")

    y -= 0.3 * inch
    c.setFont("Helvetica", 11)
    for n, entry in cm.entries():
        y = _pdf_draw_wrapped(c, 1 * inch, y, f"[{n}] {entry}")
        y -= 2

    c.showPage()
    c.save()
    buf.seek(0)
    return buf


# =========================
# VISUALS: SKY + CONSTELLATIONS + PLANETS
# =========================
def altaz_to_xy(alt_deg, az_deg, w, h):
    x = int((az_deg / 360.0) * w)
    alt_clamped = max(0.0, min(90.0, float(alt_deg)))
    y = int(h - (alt_clamped / 90.0) * (h * 0.85) - (h * 0.05))
    return x, y


def mag_to_radius(mag: float):
    # brighter star => larger radius; clamp for safety
    # typical mag range for bright catalog: [-1.5 .. 6]
    m = float(mag)
    # invert: smaller mag = brighter
    r = 3.2 - (m * 0.35)
    return max(1.0, min(4.5, r))


def mag_to_brightness(mag: float):
    # smaller mag => brighter; clamp to [120..255]
    m = float(mag)
    b = 255 - int((m + 1.5) * 18)
    return max(120, min(255, b))


def draw_glow(draw, x, y, base_radius=4, glow_radius=14, color=(255, 215, 0)):
    for r in range(glow_radius, base_radius, -2):
        alpha = int(18 * (r / glow_radius))
        c = tuple(int(color[i] * (alpha / 255)) for i in range(3))
        draw.ellipse((x - r, y - r, x + r, y + r), fill=c)
    draw.ellipse((x - base_radius, y - base_radius, x + base_radius, y + base_radius), fill=color)


def generate_sky_image(
    obs_time: Time,
    location: EarthLocation,
    visible_planets,
    visible_stars,
    show_constellations=True,
    label_stars=True,
    max_render_stars=300,
    caption_text="Simulated night sky based on your birth location",
    badge_text="",
):
    W, H = 1200, 680
    img = Image.new("RGB", (W, H), (0, 0, 0))
    draw = ImageDraw.Draw(img)

    # Aesthetic starfield background
    seed = int(obs_time.jd * 10) % 10_000_000
    rng = random.Random(seed)
    for _ in range(1400):
        x = rng.randint(0, W - 1)
        y = rng.randint(0, H - 1)
        b = int(55 + (rng.random() ** 0.45) * 200)
        r = 1 if rng.random() < 0.95 else 2
        draw.ellipse((x - r, y - r, x + r, y + r), fill=(b, b, b))

    # Render computed stars (LOD)
    if visible_stars:
        # render highest-altitude first
        stars = visible_stars[:max_render_stars]
        for s in stars:
            x, y = altaz_to_xy(s["alt"], s["az"], W, H)
            r = mag_to_radius(s.get("mag", 2.0))
            b = mag_to_brightness(s.get("mag", 2.0))
            draw.ellipse((x - r, y - r, x + r, y + r), fill=(b, b, b))

        if label_stars:
            label_color = (160, 190, 255)
            for s in stars[:7]:
                x, y = altaz_to_xy(s["alt"], s["az"], W, H)
                draw.text((x + 8, y - 10), s["name"], fill=label_color)

    # Constellations (embedded line model)
    if show_constellations:
        line_color = (90, 140, 255)
        label_color = (150, 180, 255)

        for cname, segs in CONSTELLATION_LINES.items():
            for a, b in segs:
                if a not in BRIGHT_STARS or b not in BRIGHT_STARS:
                    continue

                a_ra, a_dec, _ = BRIGHT_STARS[a]
                b_ra, b_dec, _ = BRIGHT_STARS[b]

                a_alt, a_az = project_star_to_altaz(a_ra, a_dec, obs_time, location)
                b_alt, b_az = project_star_to_altaz(b_ra, b_dec, obs_time, location)

                if a_alt <= 0 and b_alt <= 0:
                    continue

                ax, ay = altaz_to_xy(a_alt, a_az, W, H)
                bx, by = altaz_to_xy(b_alt, b_az, W, H)
                draw.line((ax, ay, bx, by), fill=line_color, width=1)

            # label constellation near first segment if visible
            if segs:
                sname = segs[0][0]
                if sname in BRIGHT_STARS:
                    ra, dec, _ = BRIGHT_STARS[sname]
                    alt, az = project_star_to_altaz(ra, dec, obs_time, location)
                    if alt > 0:
                        x, y = altaz_to_xy(alt, az, W, H)
                        draw.text((x + 8, y - 12), cname, fill=label_color)

    # Planet markers
    for p in visible_planets:
        x, y = altaz_to_xy(p["alt"], p["az"], W, H)
        draw_glow(draw, x, y, base_radius=4, glow_radius=18, color=(255, 215, 0))
        draw.text((x + 10, y - 10), p["name"], fill=(255, 255, 255))

    # Bottom overlay caption
    overlay = Image.new("RGBA", img.size, (0, 0, 0, 0))
    o = ImageDraw.Draw(overlay)
    o.rectangle((0, H - 92, W, H), fill=(0, 0, 0, 175))

    if badge_text:
        o.text((24, H - 82), badge_text, fill=(255, 220, 140))
    o.text((24, H - 56), caption_text, fill=(255, 255, 255))

    img = Image.alpha_composite(img.convert("RGBA"), overlay).convert("RGB")
    return img


def create_share_card_ig(base_img: Image.Image, birth_date: date, city: str, moment_label: str, tz_name: str,
                         planets_text: str, moon_text: str):
    img = base_img.convert("RGB")
    if img.width != img.height:
        side = min(img.width, img.height)
        left = (img.width - side) // 2
        top = (img.height - side) // 2
        img = img.crop((left, top, left + side, top + side))
    img = img.resize((1080, 1080))

    draw = ImageDraw.Draw(img)
    draw.rectangle((0, 760, 1080, 1080), fill=(0, 0, 0))

    try:
        font_big = ImageFont.truetype("DejaVuSans.ttf", 44)
        font_mid = ImageFont.truetype("DejaVuSans.ttf", 30)
        font_small = ImageFont.truetype("DejaVuSans.ttf", 24)
    except Exception:
        font_big = font_mid = font_small = None

    draw.text((40, 790), "NASA Birthday Sky", fill="white", font=font_big)
    draw.text((40, 850), f"{birth_date.isoformat()} • {city}", fill="white", font=font_mid)
    draw.text((40, 895), f"{moment_label} • {tz_name}", fill=(200, 200, 200), font=font_small)
    draw.text((40, 935), f"Planets: {planets_text}", fill=(200, 200, 200), font=font_small)
    draw.text((40, 970), f"Moon: {moon_text}", fill=(200, 200, 200), font=font_small)
    draw.text((40, 1030), "Powered by NASA Open APIs • Astropy", fill=(150, 150, 150), font=font_small)

    buf = io.BytesIO()
    img.save(buf, format="PNG")
    buf.seek(0)
    return buf


def build_share_url(birth_date: date, city: str, mode_key: str, time_hhmm: str):
    q_city = city.replace(" ", "%20")
    qs = f"?date={birth_date.isoformat()}&city={q_city}&mode={mode_key}&time={time_hhmm}&autorun=1"
    if APP_BASE_URL:
        return APP_BASE_URL.rstrip("/") + "/" + qs
    return qs


# =========================
# EDUCATOR MODE (UPGRADED)
# =========================
def educator_cards_intro(tz_name: str):
    st.markdown(
        f"""
<div class="card">
  <span class="badge">Educator Mode</span>
  <h4>What students learn</h4>
  <p><b>Coordinate systems</b> (RA/Dec → Alt/Az), <b>time standards</b> (local time → UTC), and the
  <b>visibility rule</b> (altitude &gt; 0°).</p>
  <p style="margin-top:10px;">Your location’s timezone is: <b>{tz_name}</b></p>
</div>
""",
        unsafe_allow_html=True,
    )


def render_educator_cards(tz_name: str):
    st.markdown("### 🧑🏽‍🏫 Guided Lesson Cards")

    educator_cards_intro(tz_name)

    c1, c2 = st.columns(2)
    with c1:
        st.markdown(
            """
<div class="card" style="margin-top:12px;">
  <span class="badge">Activity 1</span>
  <h4>Horizon Test (Altitude)</h4>
  <p><b>Altitude &gt; 0°</b> → above the horizon (visible). <b>Altitude ≤ 0°</b> → below the horizon.</p>
  <p style="margin-top:10px;"><b>Prompt:</b> Which planet has the highest altitude? What would you expect to see first after sunset?</p>
</div>
""",
            unsafe_allow_html=True,
        )

        st.markdown(
            """
<div class="card" style="margin-top:12px;">
  <span class="badge">Activity 3</span>
  <h4>Direction Lab (Azimuth)</h4>
  <p>Azimuth is compass direction in degrees:</p>
  <p style="margin-top:6px;">0° = North • 90° = East • 180° = South • 270° = West</p>
  <p style="margin-top:10px;"><b>Prompt:</b> If azimuth is 285°, what direction is that? (Hint: between West and North)</p>
</div>
""",
            unsafe_allow_html=True,
        )

    with c2:
        st.markdown(
            """
<div class="card" style="margin-top:12px;">
  <span class="badge">Activity 2</span>
  <h4>Moon Phase + Illumination</h4>
  <p>We estimate illumination using the Sun–Moon separation angle.</p>
  <p style="margin-top:10px;"><b>Prompt:</b> As illumination increases, how does the Moon’s appearance change? When do we call it “Full”?</p>
</div>
""",
            unsafe_allow_html=True,
        )

        st.markdown(
            f"""
<div class="card" style="margin-top:12px;">
  <span class="badge">Challenge</span>
  <h4>Timezone Reality Check</h4>
  <p>Your city uses <b>{tz_name}</b>. Re-run the same date with <b>Noon</b> vs <b>Midnight</b>.</p>
  <p style="margin-top:10px;"><b>Prompt:</b> Do altitudes change? Why does UTC matter in astronomy?</p>
</div>
""",
            unsafe_allow_html=True,
        )

    st.markdown("### ✅ Mini-Quiz (2 minutes)")
    q1 = st.radio("1) What does altitude tell you?", ["Brightness", "Height above the horizon", "Distance to the star"], index=1)
    q2 = st.radio("2) What is the global time standard used in astronomy?", ["Local time", "UTC", "Daylight Saving Time"], index=1)
    q3 = st.radio("3) Azimuth 90° points toward:", ["North", "East", "South"], index=1)

    score = int(q1 == "Height above the horizon") + int(q2 == "UTC") + int(q3 == "East")
    st.markdown(f"""
<div class="kpi">
  <div class="pill">Quiz score: <b>{score}/3</b></div>
  <div class="pill">Tip: Re-run with different cities (same date/time) and compare the sky.</div>
</div>
""", unsafe_allow_html=True)

    with st.expander("📚 Glossary (quick)", expanded=False):
        st.write(
            "- **RA/Dec**: celestial coordinates fixed on the sky (like longitude/latitude but on the celestial sphere)\n"
            "- **Alt/Az**: local sky coordinates relative to your location (what you actually see)\n"
            "- **UTC**: universal time standard used for precise astronomy\n"
            "- **Magnitude (mag)**: brightness scale (smaller/negative = brighter)\n"
        )


# =========================
# RUN ACTION
# =========================
run_clicked = st.button("🚀 Run my birthday")
should_run = run_clicked or default_autorun

if not should_run:
    st.caption("Enter your birthday + city, then press **Run my birthday** 🚀")
    st.stop()

# Validate city
if not birth_city.strip():
    st.warning('Please enter a birth city (example: "Kinshasa, Congo" or "Salt Lake City, UT").')
    st.stop()

# Geocode
with st.spinner("Finding your city coordinates..."):
    geo = geocode_city(birth_city)

if not geo:
    st.error('Couldn’t find that city. Try adding country/state (e.g., "Kinshasa, Congo").')
    st.stop()

lat, lon, resolved_address = geo
location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg)

# Moment selection
if moment.startswith("Midnight"):
    mode_key = "midnight"
    moment_label = "Midnight"
    birth_time_final = dtime(0, 0)
elif moment.startswith("Noon"):
    mode_key = "noon"
    moment_label = "Noon"
    birth_time_final = dtime(12, 0)
else:
    mode_key = "exact"
    moment_label = "Exact time"
    birth_time_final = birth_time

tz_name, utc_dt = local_dt_to_utc(birth_date, birth_time_final, lat, lon)
if not tz_name or not utc_dt:
    st.error("Could not determine timezone for this location.")
    st.stop()

obs_time = Time(utc_dt)

# Planets + Moon
with st.spinner("Calculating planets + moon (timezone-correct)..."):
    visible_planets = compute_visible_planets(location, obs_time)
    moon_illum, moon_phase_name, waxing = compute_moon_phase(obs_time)

moon_pct = int(round(moon_illum * 100))
moon_text = f"{moon_phase_name} • {moon_pct}% lit"

# Stars
obs_time_iso = utc_dt.isoformat()

if catalog_mode.startswith("CSV"):
    if DEFAULT_CSV_PATH.exists():
        with st.spinner("Computing stars from CSV catalog..."):
            visible_stars = compute_visible_stars_csv(
                obs_time_iso, lat, lon,
                csv_path_str=DEFAULT_CSV_PATH.as_posix(),
                max_rows=max_compute_stars,
            )
        star_source_used = f"CSV catalog: {DEFAULT_CSV_PATH.as_posix()} (local file), transformed to Alt/Az via Astropy."
        if not visible_stars:
            visible_stars = compute_visible_stars_embedded(obs_time_iso, lat, lon)
            star_source_used += " (CSV empty/invalid; fell back to embedded subset.)"
        catalog_mode_used = "CSV (local)"
    else:
        visible_stars = compute_visible_stars_embedded(obs_time_iso, lat, lon)
        star_source_used = "Embedded bright-star subset (CSV missing; fallback)."
        catalog_mode_used = "Embedded (fallback)"
else:
    visible_stars = compute_visible_stars_embedded(obs_time_iso, lat, lon)
    star_source_used = "Embedded bright-star subset."
    catalog_mode_used = "Embedded"

# KPIs
st.markdown(
    f"""
<div class="kpi">
  <div class="pill">Timezone: <b>{tz_name}</b></div>
  <div class="pill">UTC used: <b>{utc_dt.isoformat()}</b></div>
  <div class="pill">Stars visible: <b>{len(visible_stars)}</b></div>
</div>
""",
    unsafe_allow_html=True,
)

# Planets section
st.subheader("🪐 Planets visible above the horizon")
if visible_planets:
    ph = st.empty()
    lines = []
    for p in visible_planets:
        lines.append(f"**{p['name']}** — altitude **{p['alt']:.1f}°**, azimuth **{p['az']:.1f}°**")
        ph.markdown("\n\n".join(lines))
        time.sleep(0.18)
else:
    st.write("No major planets were above the horizon at that moment (Mars/Venus/Jupiter list).")

# Moon section
st.subheader("🌕 Moon")
st.write(f"**{moon_phase_name}** — **{moon_pct}%** illuminated")

# Educator Mode
if educator_mode:
    with st.expander("🧪 Educator Mode — explanation + lesson plan", expanded=True):
        st.markdown(
            f"""
**1) Location → Coordinates**  
We convert your city into latitude/longitude using geocoding.

**2) Coordinates → Timezone → UTC**  
Your input time is local to the birth city.  
We determine the timezone (**{tz_name}**) and convert local time to **UTC** (standard for astronomy).

**3) Planet/star visibility rule**  
- **Altitude > 0°** → visible above horizon  
- **Altitude ≤ 0°** → below horizon

**4) Azimuth is direction**  
- 0° = North, 90° = East, 180° = South, 270° = West

**5) RA/Dec → Alt/Az**  
Stars are stored in a catalog using RA/Dec, then transformed to your local sky (Alt/Az) using Astropy.
"""
        )
        render_educator_cards(tz_name)

st.divider()

# NASA APOD
with st.spinner("Fetching NASA APOD (may be degraded during outages)..."):
    apod = fetch_apod(birth_date.isoformat(), NASA_API_KEY)

apod_img = None
apod_title = ""
apod_expl = ""
apod_ok = False
apod_status_text = ""
visual_for_card = None

if isinstance(apod, dict):
    apod_title = str(apod.get("title", "")).strip()
    apod_expl = str(apod.get("explanation", "")).strip()
    apod_url = str(apod.get("url", "")).strip()
    media_type = str(apod.get("media_type", "image")).lower()

    if media_type == "image" and apod_url:
        try:
            rr = requests.get(apod_url, timeout=(4, 10))
            rr.raise_for_status()
            apod_img = Image.open(io.BytesIO(rr.content)).convert("RGB")
            apod_ok = True
        except Exception:
            apod_ok = False
    else:
        # video or unknown
        apod_ok = False
        apod_status_text = f"APOD returned media_type='{media_type}' (not an image)."

# Visual section
if apod_ok and apod_img is not None:
    st.success("✅ NASA APOD loaded successfully.")
    st.subheader(apod_title or "NASA Astronomy Picture of the Day")
    st.image(apod_img, use_container_width=True)
    if apod_expl:
        with st.expander("APOD explanation"):
            st.write(apod_expl)
    visual_for_card = apod_img
    apod_status_text = "APOD loaded successfully (NASA APOD image used)."
    badge_text = ""
    caption_text = "NASA APOD image (official NASA Open API)"
else:
    st.warning("🛑 NASA APOD is unavailable (or not an image). Planetary calculations are still accurate.")
    badge_text = "APOD unavailable — simulated sky remains accurate"
    caption_text = "Simulated sky based on your birth location (planets + moon are real calculations)"
    sky_img = generate_sky_image(
        obs_time=obs_time,
        location=location,
        visible_planets=visible_planets,
        visible_stars=visible_stars,
        show_constellations=show_constellations,
        label_stars=label_stars,
        max_render_stars=max_render_stars,
        caption_text=caption_text,
        badge_text=badge_text,
    )
    st.image(
        sky_img,
        use_container_width=True,
        caption="Simulated night sky based on your birth location (with constellation overlays if enabled)",
    )
    visual_for_card = sky_img
    if not apod_status_text:
        apod_status_text = "APOD unavailable at runtime; app displayed simulated sky visualization while preserving planet/moon calculations."

# Share URL + downloads
planets_text = ", ".join([p["name"] for p in visible_planets]) if visible_planets else "None"
time_hhmm = f"{birth_time_final.hour:02d}:{birth_time_final.minute:02d}"
share_url = build_share_url(birth_date, birth_city, mode_key, time_hhmm)

st.subheader("🔗 Share this result")
st.code(share_url, language="text")
if not APP_BASE_URL:
    st.caption("Tip: add APP_BASE_URL in Streamlit Secrets for a full clickable link.")

st.subheader("📸 Download IG Share Card")
card_buf = create_share_card_ig(
    base_img=visual_for_card,
    birth_date=birth_date,
    city=birth_city,
    moment_label=moment_label,
    tz_name=tz_name,
    planets_text=planets_text,
    moon_text=moon_text,
)
st.download_button(
    "Download PNG",
    data=card_buf,
    file_name="nasa_birthday_sky.png",
    mime="image/png",
    use_container_width=True,
)

st.subheader("📄 Export Science Report (PDF)")
pdf_buf = generate_science_report_buffer(
    birth_date_str=birth_date.isoformat(),
    city=birth_city,
    lat=lat,
    lon=lon,
    timezone=tz_name,
    utc_dt_str=utc_dt.isoformat(),
    moment_label=moment_label,
    planets=visible_planets,
    moon_text=moon_text,
    repro_link=share_url,
    apod_status_text=apod_status_text,
    apod_title=apod_title,
    apod_expl=apod_expl,
    stars_used=visible_stars,
    star_source_text=f"{STAR_CATALOG_SOURCE} | {star_source_used}",
    catalog_mode_used=catalog_mode_used,
    include_apod_excerpt=include_apod_excerpt,
)
st.download_button(
    "Download PDF Report",
    data=pdf_buf,
    file_name="nasa_birthday_sky_science_report.pdf",
    mime="application/pdf",
    use_container_width=True,
)

st.caption(f"📍 Location used: {resolved_address} (lat {lat:.3f}, lon {lon:.3f}) • Time used (UTC): {utc_dt.isoformat()}")

# Persist params so refresh keeps result
try:
    st.query_params.update({
        "date": birth_date.isoformat(),
        "city": birth_city,
        "mode": mode_key,
        "time": time_hhmm,
        "autorun": "1",
    })
except Exception:
    st.experimental_set_query_params(
        date=birth_date.isoformat(),
        city=birth_city,
        mode=mode_key,
        time=time_hhmm,
        autorun="1",
    )

