import io
import csv
import math
import time
import random
import textwrap
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
.block-container { padding-top: 2rem; padding-bottom: 3rem; max-width: 900px; }
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
.card p { margin: 0; color: rgba(255,255,255,0.82); }
.badge {
  display:inline-block;
  padding: 4px 10px;
  border-radius:999px;
  background: rgba(255,255,255,0.08);
  border: 1px solid rgba(255,255,255,0.14);
  font-size: 12px;
  margin-right: 8px;
}
</style>
""",
    unsafe_allow_html=True,
)

NASA_API_KEY = st.secrets.get("NASA_API_KEY")
APP_BASE_URL = st.secrets.get("APP_BASE_URL", "").strip()

# Optional: star catalog citation source text (you can set this in Secrets)
# Example:
# STAR_CATALOG_SOURCE = "Bright stars subset (RA/Dec) compiled for demo use; recommend Hipparcos/Tycho-2 for production."
STAR_CATALOG_SOURCE = st.secrets.get(
    "STAR_CATALOG_SOURCE",
    "Embedded bright-star subset (RA/Dec). For production-grade catalogs, use Hipparcos / Tycho-2 or Gaia DR3."
)

if not NASA_API_KEY:
    st.error('NASA API key missing. Add it in Streamlit → App settings → Secrets:\n\nNASA_API_KEY = "YOUR_KEY"')
    st.stop()

st.title("🌌 NASA Birthday Sky")
st.caption("Personal astronomy using NASA open data + real sky calculations (timezone-correct).")


# =========================
# LIGHTWEIGHT "REAL STAR CATALOG" (BRIGHT STAR SUBSET)
# =========================
# NOTE: Compact embedded subset with real RA/Dec degrees.
# Designed for reliability on Streamlit Cloud.
BRIGHT_STARS = {
    # Orion region
    "Betelgeuse": (88.7929, 7.4071),
    "Rigel": (78.6345, -8.2016),
    "Bellatrix": (81.2828, 6.3497),
    "Saiph": (86.9391, -9.6696),
    "Alnitak": (85.1897, -1.9426),
    "Alnilam": (84.0534, -1.2019),
    "Mintaka": (83.0017, -0.2991),

    # Ursa Major (Big Dipper)
    "Dubhe": (165.9322, 61.7510),
    "Merak": (165.4603, 56.3824),
    "Phecda": (178.4577, 53.6948),
    "Megrez": (183.8565, 57.0326),
    "Alioth": (193.5073, 55.9598),
    "Mizar": (200.9814, 54.9254),
    "Alkaid": (206.8856, 49.3133),

    # Leo (a few)
    "Regulus": (152.0929, 11.9672),
    "Denebola": (177.2649, 14.5721),
    "Algieba": (154.9931, 19.8415),

    # Summer triangle
    "Vega": (279.2347, 38.7837),
    "Deneb": (310.3579, 45.2803),
    "Altair": (297.6958, 8.8683),

    # Others
    "Sirius": (101.2875, -16.7161),
    "Canopus": (95.9879, -52.6957),
    "Arcturus": (213.9153, 19.1824),
    "Spica": (201.2983, -11.1613),
    "Antares": (247.3519, -26.4320),
    "Capella": (79.1723, 45.9979),
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
# STAR CATALOG LOADER (OPTIONAL CSV)
# =========================
def load_star_catalog(path="data/bright_stars.csv", max_rows=1500):
    """
    Optional CSV catalog loader.
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
                name = row.get("name", "").strip() or f"Star {i+1}"
                mag = float(row.get("mag", 0)) if row.get("mag") not in (None, "") else 0.0
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

# Optional “catalog mode” switch (still lightweight unless you add CSV)
use_csv_catalog = st.toggle("🔭 Real star catalog mode (CSV if available)", value=False)

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
    """
    Visible stars from embedded BRIGHT_STARS.
    Cached by time+location.
    """
    obs_time = Time(obs_time_iso)
    location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg)

    visible = []
    for name, (ra, dec) in BRIGHT_STARS.items():
        alt, az = project_star_to_altaz(ra, dec, obs_time, location)
        if alt > 0:
            visible.append({"name": name, "ra": ra, "dec": dec, "alt": float(alt), "az": float(az), "mag": 0.0})
    visible.sort(key=lambda s: s["alt"], reverse=True)
    return visible


@st.cache_data(show_spinner=False, ttl=60 * 60)
def compute_visible_stars_csv(obs_time_iso: str, lat: float, lon: float, max_stars=900):
    """
    Visible stars from optional CSV catalog.
    If no CSV exists, returns [].
    Cached by time+location.
    """
    catalog = load_star_catalog("data/bright_stars.csv", max_rows=1500)
    if not catalog:
        return []

    obs_time = Time(obs_time_iso)
    location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg)

    visible = []
    # Limit to keep Streamlit snappy
    for star in catalog[:max_stars]:
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
# SCIENCE REPORT (PDF EXPORT) — WITH AUTO CITATIONS
# =========================
def _pdf_wrapped_text(c, x, y, text, max_width_chars=95, line_height=14):
    for line in textwrap.wrap(text, width=max_width_chars):
        c.drawString(x, y, line)
        y -= line_height
    return y


def generate_science_report_buffer(
    birth_date,
    city,
    lat,
    lon,
    timezone,
    utc_dt,
    planets,
    moon_text,
    repro_link,
    apod_status_text,
    stars_used,
    star_source_text,
):
    buf = io.BytesIO()
    c = canvas.Canvas(buf, pagesize=LETTER)
    width, height = LETTER

    y = height - 1 * inch
    c.setFont("Helvetica-Bold", 16)
    c.drawString(1 * inch, y, "NASA Birthday Sky — Science Report")

    y -= 0.5 * inch
    c.setFont("Helvetica", 11)
    c.drawString(1 * inch, y, f"Generated: {datetime.utcnow().isoformat()} UTC")

    # Inputs
    y -= 0.55 * inch
    c.setFont("Helvetica-Bold", 12)
    c.drawString(1 * inch, y, "Inputs")

    y -= 0.3 * inch
    c.setFont("Helvetica", 11)
    c.drawString(1 * inch, y, f"Birth date: {birth_date}")
    y -= 0.25 * inch
    c.drawString(1 * inch, y, f"City: {city}")
    y -= 0.25 * inch
    c.drawString(1 * inch, y, f"Lat/Lon: {lat:.4f}, {lon:.4f}")
    y -= 0.25 * inch
    c.drawString(1 * inch, y, f"Timezone: {timezone}")
    y -= 0.25 * inch
    c.drawString(1 * inch, y, f"UTC time used: {utc_dt}")

    # Results
    y -= 0.4 * inch
    c.setFont("Helvetica-Bold", 12)
    c.drawString(1 * inch, y, "Results")

    y -= 0.3 * inch
    c.setFont("Helvetica", 11)
    c.drawString(1 * inch, y, f"Moon: {moon_text}")

    for p in planets:
        y -= 0.25 * inch
        c.drawString(1 * inch, y, f"{p['name']} — alt {p['alt']:.2f}°, az {p['az']:.2f}°")

    # Stars summary (auto)
    y -= 0.45 * inch
    c.setFont("Helvetica-Bold", 12)
    c.drawString(1 * inch, y, "Stars Used (for sky/constellation overlay)")

    y -= 0.3 * inch
    c.setFont("Helvetica", 11)
    c.drawString(1 * inch, y, f"Visible stars (sample): {min(len(stars_used), 12)} of {len(stars_used)}")

    sample = stars_used[:12]
    for s in sample:
        y -= 0.22 * inch
        c.drawString(
            1 * inch,
            y,
            f"{s['name']} — alt {s['alt']:.1f}°, az {s['az']:.1f}° (RA {s['ra']:.2f}°, Dec {s['dec']:.2f}°)",
        )

    # APOD status
    y -= 0.45 * inch
    c.setFont("Helvetica-Bold", 12)
    c.drawString(1 * inch, y, "NASA APOD Status")

    y -= 0.3 * inch
    c.setFont("Helvetica", 11)
    y = _pdf_wrapped_text(c, 1 * inch, y, apod_status_text, max_width_chars=95, line_height=14)

    # Repro
    y -= 0.35 * inch
    c.setFont("Helvetica-Bold", 12)
    c.drawString(1 * inch, y, "Reproducibility Link")

    y -= 0.3 * inch
    c.setFont("Helvetica", 11)
    y = _pdf_wrapped_text(c, 1 * inch, y, repro_link, max_width_chars=95, line_height=14)

    # Citations (auto embed)
    y -= 0.35 * inch
    c.setFont("Helvetica-Bold", 12)
    c.drawString(1 * inch, y, "Citations / Data Sources")

    y -= 0.3 * inch
    c.setFont("Helvetica", 11)
    citations = [
        "NASA Open APIs — APOD (Astronomy Picture of the Day) and related services: https://api.nasa.gov",
        "Astropy — astronomical coordinate transforms and ephemerides: https://www.astropy.org",
        f"Star catalog source: {star_source_text}",
        "Geocoding: OpenStreetMap Nominatim (Geopy client).",
        "Timezone lookup: timezonefinder + pytz.",
    ]
    for item in citations:
        y = _pdf_wrapped_text(c, 1 * inch, y, f"• {item}", max_width_chars=95, line_height=14)
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
    caption_text="Simulated night sky based on your birth location",
    badge_text="",
):
    W, H = 1080, 620
    img = Image.new("RGB", (W, H), (0, 0, 0))
    draw = ImageDraw.Draw(img)

    # Base starfield (aesthetic filler)
    seed = int(obs_time.jd * 10) % 10_000_000
    rng = random.Random(seed)
    for _ in range(1200):
        x = rng.randint(0, W - 1)
        y = rng.randint(0, H - 1)
        b = int(70 + (rng.random() ** 0.4) * 185)
        r = 1 if rng.random() < 0.95 else 2
        draw.ellipse((x - r, y - r, x + r, y + r), fill=(b, b, b))

    # Overlay real visible star points (intentional “catalog vibe”)
    # Brighter stars (higher altitude) get slightly larger
    if visible_stars:
        for s in visible_stars[:220]:  # cap for speed/clarity
            x, y = altaz_to_xy(s["alt"], s["az"], W, H)
            r = 1 if s["alt"] < 20 else 2
            b = 210 if s["alt"] > 45 else 180
            draw.ellipse((x - r, y - r, x + r, y + r), fill=(b, b, b))

        # Label a few highest-altitude stars (subtle)
        label_color = (160, 190, 255)
        for s in visible_stars[:6]:
            x, y = altaz_to_xy(s["alt"], s["az"], W, H)
            draw.text((x + 8, y - 10), s["name"], fill=label_color)

    # Constellations (from embedded set)
    if show_constellations:
        line_color = (90, 140, 255)
        label_color = (150, 180, 255)

        for cname, segs in CONSTELLATION_LINES.items():
            for a, b in segs:
                if a not in BRIGHT_STARS or b not in BRIGHT_STARS:
                    continue

                a_ra, a_dec = BRIGHT_STARS[a]
                b_ra, b_dec = BRIGHT_STARS[b]

                a_alt, a_az = project_star_to_altaz(a_ra, a_dec, obs_time, location)
                b_alt, b_az = project_star_to_altaz(b_ra, b_dec, obs_time, location)

                if a_alt <= 0 and b_alt <= 0:
                    continue

                ax, ay = altaz_to_xy(a_alt, a_az, W, H)
                bx, by = altaz_to_xy(b_alt, b_az, W, H)
                draw.line((ax, ay, bx, by), fill=line_color, width=1)

            first_seg = segs[0] if segs else None
            if first_seg:
                sname = first_seg[0]
                if sname in BRIGHT_STARS:
                    ra, dec = BRIGHT_STARS[sname]
                    alt, az = project_star_to_altaz(ra, dec, obs_time, location)
                    if alt > 0:
                        x, y = altaz_to_xy(alt, az, W, H)
                        draw.text((x + 8, y - 12), cname, fill=label_color)

    # Planet markers
    for p in visible_planets:
        x, y = altaz_to_xy(p["alt"], p["az"], W, H)
        draw_glow(draw, x, y, base_radius=4, glow_radius=16, color=(255, 215, 0))
        draw.text((x + 10, y - 10), p["name"], fill=(255, 255, 255))

    # Bottom overlay caption (designed, intentional)
    overlay = Image.new("RGBA", img.size, (0, 0, 0, 0))
    o = ImageDraw.Draw(overlay)
    o.rectangle((0, H - 86, W, H), fill=(0, 0, 0, 170))

    # “badge”
    if badge_text:
        o.text((24, H - 78), badge_text, fill=(255, 220, 140))

    o.text((24, H - 54), caption_text, fill=(255, 255, 255))
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
# EDUCATOR MODE LESSON CARDS
# =========================
def render_educator_cards(tz_name: str):
    st.markdown("### 🧑🏽‍🏫 Educator Mode — Guided Lesson Cards")
    c1, c2 = st.columns(2)
    with c1:
        st.markdown(
            """
<div class="card">
  <span class="badge">Activity 1</span>
  <h4>Horizon Test (Altitude)</h4>
  <p>Rule: <b>Altitude &gt; 0°</b> means “above the horizon” (visible).</p>
  <p>Ask: Which planet has the highest altitude? What would you expect to see first after sunset?</p>
</div>
""",
            unsafe_allow_html=True,
        )
        st.markdown(
            """
<div class="card" style="margin-top:12px;">
  <span class="badge">Activity 3</span>
  <h4>Direction Lab (Azimuth)</h4>
  <p>Azimuth is direction: 0°=N, 90°=E, 180°=S, 270°=W.</p>
  <p>Ask: If azimuth is 285°, what compass direction is that?</p>
</div>
""",
            unsafe_allow_html=True,
        )

    with c2:
        st.markdown(
            """
<div class="card">
  <span class="badge">Activity 2</span>
  <h4>Moon Phase + Illumination</h4>
  <p>We estimate illumination from Sun–Moon separation.</p>
  <p>Ask: How does the moon’s appearance change as illumination increases?</p>
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
  <p>Ask: Do planet altitudes change? Why?</p>
</div>
""",
            unsafe_allow_html=True,
        )

    st.caption("Tip: Have students compare two cities on the same date/time and explain differences using latitude + timezone.")


# =========================
# RUN ACTION
# =========================
run_clicked = st.button("🚀 Run my birthday")
should_run = run_clicked or default_autorun

if should_run:
    if not birth_city.strip():
        st.warning("Please enter a birth city (example: “Kinshasa, Congo” or “Salt Lake City, UT”).")
        st.stop()

    with st.spinner("Finding your city coordinates..."):
        geo = geocode_city(birth_city)

    if not geo:
        st.error("Couldn’t find that city. Try adding country/state (e.g., “Kinshasa, Congo”).")
        st.stop()

    lat, lon, resolved_address = geo
    location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg)

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

    # compute planets + moon
    with st.spinner("Calculating planets + moon (timezone-correct)..."):
        visible_planets = compute_visible_planets(location, obs_time)
        moon_illum, moon_phase_name, waxing = compute_moon_phase(obs_time)

    moon_pct = int(round(moon_illum * 100))
    moon_text = f"{moon_phase_name} • {moon_pct}% lit"

    # Stars: embedded or CSV
    obs_time_iso = utc_dt.isoformat()
    if use_csv_catalog:
        visible_stars = compute_visible_stars_csv(obs_time_iso, lat, lon)
        star_source_used = "CSV catalog: data/bright_stars.csv (local file), transformed to Alt/Az via Astropy."
        if not visible_stars:
            visible_stars = compute_visible_stars_embedded(obs_time_iso, lat, lon)
            star_source_used = "Embedded bright-star subset (fallback)."
    else:
        visible_stars = compute_visible_stars_embedded(obs_time_iso, lat, lon)
        star_source_used = "Embedded bright-star subset."

    # Planet “animation sweep” (nice vibe)
    st.subheader("🪐 Planets visible above the horizon")
    if visible_planets:
        ph = st.empty()
        lines = []
        for p in visible_planets:
            lines.append(f"**{p['name']}** — altitude **{p['alt']:.1f}°**, azimuth **{p['az']:.1f}°**")
            ph.markdown("\n\n".join(lines))
            time.sleep(0.25)
    else:
        st.write("No major planets were above the horizon at that moment (Mars/Venus/Jupiter list).")

    st.subheader("🌕 Moon")
    st.write(f"**{moon_phase_name}** — **{moon_pct}%** illuminated")

    # Educator mode: explain + lesson cards
    if educator_mode:
        with st.expander("🧪 Educator Mode — explanation + lesson plan", expanded=True):
            st.markdown(
                f"""
**1) Location → Coordinates**  
We convert your city into latitude/longitude using geocoding.

**2) Coordinates → Timezone → UTC**  
Your input time is local to the birth city.  
We determine the timezone (**{tz_name}**) and convert local time to **UTC** (the standard for astronomy).

**3) Planet visibility = altitude > 0°**  
Altitude tells how high an object is above the horizon:
- **Altitude > 0°** → visible above horizon  
- **Altitude ≤ 0°** → below horizon

**4) Azimuth is direction**  
Azimuth is a compass direction in degrees:
- 0° = North, 90° = East, 180° = South, 270° = West
"""
            )
            render_educator_cards(tz_name)

    # NASA APOD (optional; may be degraded)
    with st.spinner("Fetching NASA APOD (may be degraded during outages)..."):
        apod = fetch_apod(birth_date.isoformat(), NASA_API_KEY)

    apod_img = None
    apod_title = ""
    apod_expl = ""
    apod_ok = False

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

    st.divider()

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
        # Intentional, designed fallback messaging
        st.warning("🛑 NASA APOD is temporarily unavailable (NASA shows a service outage). Planetary calculations are still accurate.")
        badge_text = "APOD temporarily unavailable — simulated sky remains accurate"
        caption_text = "Simulated night sky based on your birth location (planets + moon are real calculations)"
        sky_img = generate_sky_image(
            obs_time=obs_time,
            location=location,
            visible_planets=visible_planets,
            visible_stars=visible_stars,
            show_constellations=True,
            caption_text=caption_text,
            badge_text=badge_text,
        )
        st.image(
            sky_img,
            use_container_width=True,
            caption="Simulated night sky based on your birth location (with constellation overlays)",
        )
        visual_for_card = sky_img
        apod_status_text = "APOD unavailable at runtime; app displayed simulated sky visualization while preserving planet/moon calculations."

    # Share link + download card
    planets_text = ", ".join([p["name"] for p in visible_planets]) if visible_planets else "None"
    time_hhmm = f"{birth_time_final.hour:02d}:{birth_time_final.minute:02d}"

    st.subheader("🔗 Share this result")
    share_url = build_share_url(birth_date, birth_city, mode_key, time_hhmm)
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

    # ✅ NEW: Science Report PDF export with citations + repro link
    st.subheader("📄 Export Science Report (PDF)")
    pdf_buf = generate_science_report_buffer(
        birth_date=birth_date.isoformat(),
        city=birth_city,
        lat=lat,
        lon=lon,
        timezone=tz_name,
        utc_dt=utc_dt.isoformat(),
        planets=visible_planets,
        moon_text=moon_text,
        repro_link=share_url,
        apod_status_text=apod_status_text,
        stars_used=visible_stars,
        star_source_text=f"{STAR_CATALOG_SOURCE} | {star_source_used}",
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

else:
    st.caption("Enter your birthday + city, then press **Run my birthday** 🚀")
